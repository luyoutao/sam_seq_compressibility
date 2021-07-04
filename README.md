# Summary:
Compute gzip compressibility for each read (read pair if paired-end [PE]). The compressibility is defined as the original string length divided by the compressed length, and the higher compressibility the lower complexity.

## Usage:

    sam_seq_compressibility --inFile input.bam [--outFile output.tsv] [--nameSorted] [--ncores 1] [--endType PE] [--include_softclipped] [--version|-v]

## Output:
The output has 7 columns:

    1) read ID;
    2) sequence of R1;
    3) sequence of R2 (if any);
    4) sequence of both mates (concatenated but overlapping region (if any) merged)
    5) compressibility of R1;
    6) compressibility of R2 (if any);
    7) compressibility of both mates (concatenated but overlapping region (if any) merged)

## Options:

    -i, --inFile input file, can be SAM or BAM format

    -o, --outFile output file; if omitted, write to STDOUT

        --endType       choose from SE or PE (default: 'PE')
        --nameSorted    whether the input files are ALL name sorted. If not,
                        they will be name sorted and saved to temporary files
                        with '.nameSorted.bam' suffix
        --ncores number of cores to use for name sorting (default: 1)

        --include_softclipped
                        keep the softclipped part from being trimmed
    -h, --help          print usage
    -v, --version       print version

## Example:

    ./sam_seq_compressibility -i test/input.bam 
    ReadID  SeqR1   SeqR2   SeqPair CompressibilityR1       CompressibilityR2       CompressibilityPair
    NB501328:197:HMK3KBGX7:1:11101:1498:13562       GGGGGGGCGGGGGGGGGGGGGGGGTGGGGGGGTGGGGTGGGGGGGGGGGGGGGGGG                GGGGGGGCGGGGGGGGGGGGGGGGTGGGGGGGTGGGGTGGGGGGGGGGGGGGGGGG        2.28    1       2.28
    NB501328:197:HMK3KBGX7:1:11101:1744:6598        TTTCCGTGATTTTCAGTTTTCTTT        TTTCCGTGATTTTCAGTTTTCTTT        TTTCCGTGATTTTCAGTTTTCTTT        1.3888888       1.3888888       1.3888888
    NB501328:197:HMK3KBGX7:1:11101:2468:19747       AATGCTCAGAGAGACAGATATAAGCCTTATATTTCAGCCCTCTTCTCAACTTTA  CAAAAAAACCTCCAAACTCATTCTCATAAGTCAATGACCTACATATTCAATTG   CAAAAAAACCTCCAAACTCATTCTCATAAGTCAATGACCTACATATTCAATTGAATGCTCAGAGAGACAGATATAAGCCTTATATTTCAGCCCTCTTCTCAACTTTA 1.5714285       1.5882353       2.0769231
    NB501328:197:HMK3KBGX7:1:11101:2835:13039       AGTGGGGGGGAAAATCCTGATCAAT       AGTGGGGGGGAAAATCCTGATCAAT       AGTGGGGGGGAAAATCCTGATCAAT       1.5294118       1.5294118       1.5294118
