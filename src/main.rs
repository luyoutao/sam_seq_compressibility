#' This Source Code Form is subject to the terms of the Mozilla Public
#' License, v. 2.0. If a copy of the MPL was not distributed with this
#' file, You can obtain one at http://mozilla.org/MPL/2.0/.
#'
#' Youtao Lu@Kim Lab, 2016-2020
#' Kim Laboratory Bioinformatics Infrastructure

use chrono::Local;
use env_logger::{self, Builder};
use flate2::read::GzEncoder;
use flate2::Compression;
use getopts::Options;
use log::{debug, error, info, trace, warn, LevelFilter};
use regex::Regex;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::exit;
use std::process::{Command, Stdio};
use std::str::from_utf8;

static VERSION: &str = "0.2.2";

#[derive(Debug)]
enum EndType {
    SE,
    PE,
}

#[derive(Debug)]
struct Params {
    infile: String,
    outfile: Option<String>,
    endtype: EndType,
    namesorted: bool,
    ncores: u8,
    include_softclipped: bool,
}

#[derive(Debug)]
struct Mate {
    id: String,
    flag: u16,
    chr: String,
    pos: u32,
    cigar: String,
    seq: String,
    tags: Vec<String>,
    clip_l: Option<usize>,
    clip_r: Option<usize>,
    is_r1: Option<bool>,
    is_rev: Option<bool>,
}

impl Clone for Mate {
    fn clone(&self) -> Mate {
        Mate {
            id: self.id.clone(),
            flag: self.flag,
            chr: self.chr.clone(),
            pos: self.pos,
            cigar: self.cigar.clone(),
            seq: self.seq.clone(),
            tags: self.tags.iter().map(|s| s.clone()).collect(),
            clip_l: self.clip_l,
            clip_r: self.clip_r,
            is_r1: self.is_r1,
            is_rev: self.is_rev,
        }
    }
}

impl Mate {
    fn new(
        id: String,
        flag: u16,
        chr: String,
        pos: u32,
        cigar: String,
        seq: String,
        tags: Vec<String>,
    ) -> Option<Self> {
        Some(Self {
            id: id,
            flag: flag,
            chr: chr,
            pos: pos,
            cigar: cigar,
            seq: seq,
            tags: tags,
            clip_l: None,
            clip_r: None,
            is_r1: None,
            is_rev: None,
        })
    }

    fn set_is_r1(&mut self) {
        self.is_r1 = Some(self.flag & 0x40 != 0);
    }

    fn set_is_rev(&mut self) {
        self.is_rev = Some(self.flag & 0x10 != 0);
    }
    fn parse_flag(&mut self) {
        self.set_is_r1();
        self.set_is_rev();
    }

    fn parse_cigar(&mut self) {
        let re_l = Regex::new(r"^(\d+)S").expect("Failed to create regex for ^(\\d+)S");
        let re_r = Regex::new(r"(\d+)S$").expect("Failed to create regex for (\\d+)S$");
        let cap_l = re_l.captures(&self.cigar);
        let cap_r = re_r.captures(&self.cigar);

        let clip_l = match cap_l {
            None => 0,
            Some(c) => c
                .get(1)
                .map_or("0", |m| m.as_str())
                .parse::<usize>()
                .unwrap(),
        };

        let clip_r = match cap_r {
            None => 0,
            Some(c) => c
                .get(1)
                .map_or("0", |m| m.as_str())
                .parse::<usize>()
                .unwrap(),
        };
        self.clip_l = Some(clip_l);
        self.clip_r = Some(clip_r);
    }

    fn get_seq(&self, include_softclipped: bool) -> Option<String> {
        if include_softclipped {
            return Some(self.seq.to_string());
        } else {
            let n = self.seq.chars().count();
            let substr = self
                .seq
                .chars()
                .skip(self.clip_l.unwrap())
                .take(n - self.clip_l.unwrap() - self.clip_r.unwrap())
                .collect::<String>();
            return Some(substr);
        }
    }
}

#[derive(Debug)]
struct ReadPair {
    id: String,
    pair: [Option<Mate>; 2],
    is_discordant: Option<bool>,
    overlap: Option<usize>,
}

impl ReadPair {
    fn new(r1: Option<Mate>, r2: Option<Mate>) -> Option<Self> {
        if r1.as_ref().is_none() && r2.as_ref().is_none() {
            panic!("Neither mate exists!");
        }
        let id: String = if r1.as_ref().is_none() {
            String::from(&r2.as_ref().unwrap().id)
        } else {
            String::from(&r1.as_ref().unwrap().id)
        };
        Some(Self {
            id: id,
            pair: [r1, r2],
            is_discordant: None,
            overlap: None,
        })
    }

    fn set_is_discordant(&mut self) {
        let r1 = self.pair[0].as_ref();
        let r2 = self.pair[1].as_ref();
        if r1.is_none() || r2.is_none() {
            self.is_discordant = Some(false);
            return;
        }
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();
        if r1.chr != r2.chr {
            self.is_discordant = Some(true);
            return;
        }

        // [------]-----R1------X1[---->]
        // [<-----]-----R2----X2[-------]
        // if X1 > X2, then is discordant
        if !r1.is_rev.unwrap() {
            // if R1 is on the plus strand of the genome
            if r1.pos as usize + r1.seq.chars().count() - r1.clip_l.unwrap() - r1.clip_r.unwrap()
                > r2.pos as usize + r2.seq.chars().count() - r2.clip_l.unwrap() - r2.clip_r.unwrap()
            {
                self.is_discordant = Some(true);
            } else {
                self.is_discordant = Some(false);
            }
        } else {
            if r2.pos as usize + r2.seq.chars().count() - r2.clip_l.unwrap() - r2.clip_r.unwrap()
                > r1.pos as usize + r1.seq.chars().count() - r1.clip_l.unwrap() - r1.clip_r.unwrap()
            {
                self.is_discordant = Some(true);
            } else {
                self.is_discordant = Some(false);
            }
        }
    }

    fn cal_overlap(&mut self) {
        if self.is_discordant.unwrap() {
            self.overlap = Some(0);
            return;
        }
        let r1 = self.pair[0].as_ref();
        let r2 = self.pair[1].as_ref();
        if r1.is_none() || r2.is_none() {
            self.overlap = Some(0);
            return;
        }
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();

        if !r1.is_rev.unwrap() {
            let x = r1.pos as i64 + r1.seq.chars().count() as i64
                - r1.clip_l.unwrap() as i64
                - r1.clip_r.unwrap() as i64
                - r2.pos as i64;
            if x > 0 {
                self.overlap = Some(x as usize);
            } else {
                self.overlap = Some(0);
            }
        } else {
            let x = r2.pos as i64 + r2.seq.chars().count() as i64
                - r2.clip_l.unwrap() as i64
                - r2.clip_r.unwrap() as i64
                - r1.pos as i64;
            if x > 0 {
                self.overlap = Some(x as usize);
            } else {
                self.overlap = Some(0);
            }
        }
    }

    fn get_seq(&self, include_softclipped: bool) -> Option<String> {
        let r1 = self.pair[0].as_ref();
        let r2 = self.pair[1].as_ref();
        if r1.is_none() {
            return r2.unwrap().get_seq(include_softclipped);
        }
        if r2.is_none() {
            return r1.unwrap().get_seq(include_softclipped);
        }
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();
        let overlap = self.overlap.unwrap();
        if overlap == 0 {
            // if no overlap, we simply concatenate the two sequences: keep the sequence that maps to the plus strand of the genome first.
            let mut s1 = r1.get_seq(include_softclipped).unwrap().to_string();
            let mut s2 = r2.get_seq(include_softclipped).unwrap().to_string();
            if !r1.is_rev.unwrap() {
                s1.push_str(&s2);
                return Some(s1);
            } else {
                s2.push_str(&s1);
                return Some(s2);
            }
        } else {
            // if there is overlap, we pad the leftmost softclipped end from the read that maps to the plus stand of the genome, then the rightmost softclipped end from the read that maps to the minus stand of the genome, to the concatenated sequence
            let mut s1 = r1.get_seq(false).unwrap().to_string();
            let mut s2 = r2.get_seq(false).unwrap().to_string();
            let seq1 = &r1.seq;
            let seq2 = &r2.seq;
            if !r1.is_rev.unwrap() {
                s2 = s2
                    .chars()
                    .skip(r2.clip_l.unwrap() + overlap)
                    .collect::<String>();
                s1.push_str(&s2);
                if include_softclipped {
                    let mut pad1: String = seq1[0..r1.clip_l.unwrap()].to_string();
                    let pad2: String = seq2[(seq2.len() - r2.clip_r.unwrap())..].to_string();
                    pad1.push_str(&s1);
                    pad1.push_str(&pad2);
                    return Some(pad1);
                } else {
                    return Some(s1);
                }
            } else {
                s1 = s1
                    .chars()
                    .skip(r1.clip_l.unwrap() + overlap)
                    .collect::<String>();
                s2.push_str(&s1);
                if include_softclipped {
                    let mut pad2: String = seq2[0..r2.clip_l.unwrap()].to_string();
                    let pad1: String = seq1[(seq1.len() - r1.clip_r.unwrap())..].to_string();
                    pad2.push_str(&s2);
                    pad2.push_str(&pad1);
                    return Some(pad2);
                } else {
                    return Some(s2);
                }
            }
        }
    }
}

struct SamReader {
    infile: String,
    infh: Option<BufReader<std::process::ChildStdout>>,
    outfile: Option<String>,
    outfh: Option<Box<dyn Write>>,
    endtype: EndType,
    namesorted: bool,
    ncores: u8,
    include_softclipped: bool,
}

impl SamReader {
    fn new(params: Params) -> Self {
        Self {
            infile: params.infile,
            infh: None,
            outfile: params.outfile,
            outfh: None,
            endtype: params.endtype,
            namesorted: params.namesorted,
            ncores: params.ncores,
            include_softclipped: params.include_softclipped,
        }
    }

    fn init(&mut self) {
        if !self.namesorted {
            let t = self.infile.split(".").collect::<Vec<&str>>();
            let mut outfile = t.split_last().unwrap().1.join(".");
            outfile.push_str(".nameSorted.bam");
            namesort(&self.infile, &outfile, self.ncores);
            self.infile = outfile;
        }
        let child = match Command::new("samtools")
            .args(&["view", "-F", "0x4", "-F", "0x100", &self.infile])
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn()
        {
            Ok(r) => r,
            Err(e) => panic!("Cannot to read pipe {}!", e),
        };
        self.infh = Some(BufReader::new(child.stdout.unwrap()));
        match &self.outfile {
            Some(f) => self.outfh = Some(Box::new(File::create(f).unwrap()) as Box<dyn Write>),
            None => self.outfh = Some(Box::new(io::stdout()) as Box<dyn Write>),
        };
    }

    fn iter(&mut self) {
        let mut infh = self.infh.as_mut().expect("Failed to access infh!");
        let mut outfh = self.outfh.as_mut().expect("Failed to access outfh");
        let include_softclipped = self.include_softclipped;
        let mut prev_read: Option<Mate> = None;
        let mut read_pair: Option<ReadPair> = None;
        outfh.write(format!("ReadID\tSeqR1\tSeqR2\tSeqPair\tCompressibilityR1\tCompressibilityR2\tCompressibilityPair\n").as_bytes());
        for line in infh.lines() {
            let l = line.unwrap();
            let fields: Vec<String> = l.split("\t").map(|s| s.to_string()).collect();
            let id = fields[0].clone();
            let flag = fields[1].parse::<u16>().unwrap();
            let chr = fields[2].clone();
            let pos = fields[3].parse::<u32>().unwrap();
            let cigar = fields[5].clone();
            let seq = fields[9].clone();
            let tags: Vec<String> = fields[11..].to_vec();
            let mut curr = Mate::new(id, flag, chr, pos, cigar, seq, tags)
                .expect("Failed to instantialize read!");

            curr.parse_flag();
            curr.parse_cigar();
            trace!("{:?}", curr);
            match prev_read {
                Some(prev) => {
                    if prev.id == curr.id {
                        if curr.is_r1.unwrap() {
                            read_pair = ReadPair::new(Some(curr), Some(prev));
                        } else {
                            read_pair = ReadPair::new(Some(prev), Some(curr));
                        }
                        prev_read = None;
                    } else {
                        if prev.is_r1.unwrap() {
                            read_pair = ReadPair::new(Some(prev), None);
                        } else {
                            read_pair = ReadPair::new(None, Some(prev));
                        }
                        prev_read = Some(curr);
                    }
                    let mut rp = read_pair.unwrap();
                    rp.set_is_discordant();
                    rp.cal_overlap();
                    trace!("{:?}", rp);
                    let id = &rp.id;
                    let seq1 = if rp.pair[0].is_none() {
                        String::from("")
                    } else {
                        (&rp).pair[0]
                            .as_ref()
                            .unwrap()
                            .get_seq(include_softclipped)
                            .unwrap()
                    };
                    let seq2 = if rp.pair[1].is_none() {
                        String::from("")
                    } else {
                        (&rp).pair[1]
                            .as_ref()
                            .unwrap()
                            .get_seq(include_softclipped)
                            .unwrap()
                    };
                    let seqs = rp.get_seq(include_softclipped).unwrap();
                    let len1 = &seq1.as_bytes().len();
                    let len2 = &seq2.as_bytes().len();
                    let lens = &seqs.as_bytes().len();
                    outfh.write(format!("{}\t{}\t{}\t{}", &id, &seq1, &seq2, &seqs).as_bytes());
                    // ""   ->    [1f, 8b, 08, 00, 00, 00, 00, 00, 02, ff, 03, 00] 12-byte overhead
                    let mut _buf = [0u8; 1024];
                    let mut gz = GzEncoder::new(seq1.as_bytes(), Compression::best());
                    let gzlen1 = match gz.read(&mut _buf) {
                        Ok(x) => x,
                        Err(_e) => 0,
                    };
                    gz = GzEncoder::new(seq2.as_bytes(), Compression::best());
                    let gzlen2 = match gz.read(&mut _buf) {
                        Ok(x) => x,
                        Err(_e) => 0,
                    };
                    gz = GzEncoder::new(seqs.as_bytes(), Compression::best());
                    let gzlens = match gz.read(&mut _buf) {
                        Ok(x) => x,
                        Err(_e) => 0,
                    };
                    let ratio1: f32 = (len1 + 1) as f32 / (gzlen1 - 11) as f32;
                    let ratio2: f32 = (len2 + 1) as f32 / (gzlen2 - 11) as f32;
                    let ratio: f32 = (lens + 1) as f32 / (gzlens - 11) as f32;
                    outfh.write(format!("\t{}\t{}\t{}\n", ratio1, ratio2, ratio).as_bytes());
                }
                None => {
                    prev_read = Some(curr);
                }
            }
        }
        match prev_read {
            None => {}
            Some(prev) => {
                if prev.is_r1.unwrap() {
                    read_pair = ReadPair::new(Some(prev), None);
                } else {
                    read_pair = ReadPair::new(None, Some(prev));
                }
                let mut rp = read_pair.unwrap();
                rp.set_is_discordant();
                rp.cal_overlap();
                trace!("{:?}", rp);
                let id = &rp.id;
                let seq1 = if rp.pair[0].is_none() {
                    String::from("")
                } else {
                    (&rp).pair[0]
                        .as_ref()
                        .unwrap()
                        .get_seq(include_softclipped)
                        .unwrap()
                };
                let seq2 = if rp.pair[1].is_none() {
                    String::from("")
                } else {
                    (&rp).pair[1]
                        .as_ref()
                        .unwrap()
                        .get_seq(include_softclipped)
                        .unwrap()
                };
                let seqs = rp.get_seq(include_softclipped).unwrap();
                let len1 = &seq1.as_bytes().len();
                let len2 = &seq2.as_bytes().len();
                let lens = &seqs.as_bytes().len();
                outfh.write(format!("{}\t{}\t{}\t{}\n", &id, &seq1, &seq2, &seqs).as_bytes());
                let mut _buf = [0u8; 1024];
                let mut gz = GzEncoder::new(seq1.as_bytes(), Compression::best());
                let gzlen1 = match gz.read(&mut _buf) {
                    Ok(x) => x,
                    Err(_e) => 0,
                };
                gz = GzEncoder::new(seq2.as_bytes(), Compression::best());
                let gzlen2 = match gz.read(&mut _buf) {
                    Ok(x) => x,
                    Err(_e) => 0,
                };
                gz = GzEncoder::new(seqs.as_bytes(), Compression::best());
                let gzlens = match gz.read(&mut _buf) {
                    Ok(x) => x,
                    Err(_e) => 0,
                };
                let ratio1: f32 = (len1 + 1) as f32 / (gzlen1 - 11) as f32;
                let ratio2: f32 = (len2 + 1) as f32 / (gzlen2 - 11) as f32;
                let ratio: f32 = (lens + 1) as f32 / (gzlens - 11) as f32;
                outfh.write(format!("\t{}\t{}\t{}\n", ratio1, ratio2, ratio).as_bytes());
            }
        }
    }
}

fn init_logger() {
    Builder::new()
        .format(|buf, record| {
            writeln!(
                buf,
                "[{} {}] {}",
                Local::now().format("%Y-%m-%d %H:%M:%S%.3f %z"),
                record.level(),
                record.args()
            )
        })
        .filter(None, LevelFilter::Trace)
        .init();
}

fn usage(arg0: &str, opts: Options) {
    let s = format!("\
Summary:
    Computes gzip compressibility for each read (read pair if paired-end [PE]). The compressibility is defined as the original string length divided by the compressed length, and the higher compressibility the lower complexity. 

Usage:
    {} --inFile input.bam --outFile output.tsv [--nameSorted] [--ncores 1] [--endType PE] [--include_softclipped] [--version|-v]

Output:
    The output has 7 columns:
        *) read ID;
        *) sequence of R1;
        *) sequence of R2 (if any);
        *) sequence of both mates (concatenated but overlapping region (if any) merged)
        *) compressibility of R1;
        *) compressibility of R2 (if any);
        *) compressibility of both mates (concatenated but overlapping region (if any) merged)", 
arg0);
    eprintln!("{}", opts.usage(&s));
}

fn proc_args(args: &Vec<String>, mut opts: Options) -> Params {
    opts.optopt("i", "inFile", "", "input file, can be SAM or BAM format");
    opts.optopt(
        "o",
        "outFile",
        "",
        "output file; if omitted, write to STDOUT",
    );
    opts.optopt("", "endType", "choose from SE or PE (default: 'PE')", "");
    opts.optflag("", "nameSorted", "whether the input files are ALL name sorted. If not, they will be name sorted and saved to temporary files with '.nameSorted.bam' suffix");
    opts.optopt(
        "",
        "ncores",
        "",
        "number of cores to use for name sorting (default: 1)",
    );
    opts.optflag(
        "",
        "include_softclipped",
        "keep the softclipped part from being trimmed",
    );
    opts.optflag("h", "help", "print usage");
    opts.optflag("v", "version", "print version");
    let matches = opts.parse(&args[1..]).unwrap();
    if matches.opt_present("h") {
        usage(&args[0], opts);
        exit(0);
    }
    if matches.opt_present("v") {
        println!("{} v{}", &args[0], &VERSION);
        exit(0);
    }
    let infile = match matches.opt_str("inFile") {
        Some(f) => match Path::new(&f).exists() {
            true => match &*(f
                .split('.')
                .last()
                .expect("Faied to find the file extension!")
                .to_lowercase())
            {
                "sam" | "bam" => f,
                _ => panic!("{} does not seem to be a SAM or BAM!", f),
            },
            false => panic!("{} does not exist!", f),
        },
        None => panic!("--inFile is empty!"),
    };
    let outfile = matches.opt_str("outFile");

    let endtype = match matches.opt_get_default("endType", "PE".to_string()) {
        Ok(s) => match &*s {
            "PE" => EndType::PE,
            "SE" => EndType::SE,
            _ => panic!("--endType can only be 'PE' or 'SE'!"),
        },
        Err(e) => panic!("{} invalid --endType!", e),
    };
    let namesorted = matches.opt_present("nameSorted");
    let ncores = matches
        .opt_get_default("ncores", 1)
        .expect("Failed to parse --ncores!");
    let include_softclipped = matches.opt_present("include_softclipped");

    let params = Params {
        infile,
        outfile,
        endtype,
        namesorted,
        ncores,
        include_softclipped,
    };
    return params;
}

fn namesort(infile: &str, outfile: &str, ncores: u8) {
    let mut child = Command::new("samtools")
        .args(&[
            "sort",
            "-n",
            "-@",
            &ncores.to_string(),
            "-o",
            outfile,
            infile,
        ])
        .spawn()
        .expect("samtools sort failed to run!");
    let exit_status = child.wait().expect("The child process failed to return!");
    if !exit_status.success() {
        match exit_status.code() {
            None => panic!("Unknown exit code!"),
            Some(e) => panic!("samtools has nonzero exit code {}!", e),
        }
    }
}

fn main() {
    init_logger();
    let args: Vec<String> = env::args().collect();
    let params = proc_args(&args, Options::new());
    info!("{:?}", params);
    let mut samreader = SamReader::new(params);
    samreader.init();
    info!("Starting iteration...");
    samreader.iter();
    info!("All done!");
}
