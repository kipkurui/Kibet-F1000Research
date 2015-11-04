"""
Script to extract the negative bed coordinates given a bed file 
distance downstream

Usage:
    python extract_negative.py <Bed-file> <out-file> <downstream-distance>
"""
import os
import sys


def extract_negative(bed_in, bed_out, downstream_distance):
    """
    Extract bed file coordinates of the negative sequences
    """
    with open(bed_in) as bed:
        with open(bed_out, "w") as neg_fa:
            for line in bed:
                spl = line.split()
                write_details = "%s\t%i\t%i\n" % (spl[0], int(spl[1])+int(downstream_distance),
                                       int(spl[2])+int(downstream_distance))
                neg_fa.write(write_details)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print __doc__
        sys.exit(1)
    bed_in = sys.argv[1]
    fasta_out = sys.argv[2]
    dist = sys.argv[3]
    extract_negative(bed_in, fasta_out, dist)
