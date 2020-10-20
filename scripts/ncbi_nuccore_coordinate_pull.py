#/usr/bin/env python3.7

import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input gene query file, accession number in 9th field, start in 10th, stop in 11th')
args = parser.parse_args()

with open(args.input, "r") as f:
    for line in f:
        try:
            record = line.split("\t")
            t=record[2]
        except IndexError:
            break
        record = line.split("\t")
        accession = record[1]
        if int(record[8]) >= int(record[9]):
            fwdcoord = int(record[9])+1000
            revcoord = int(record[8])-1000
        else:
            fwdcoord = int(record[8])+1000
            revcoord = int(record[9].rstrip("\n"))-1000
        subprocess.check_call(['efetch -db nuccore -id %s \
            -seq_start %s \
            -seq_stop %s \
            -format fasta' %(accession, fwdcoord, revcoord)], shell=True )