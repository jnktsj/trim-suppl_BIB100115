#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Junko Tsuji

# This script trims the 3' end of each read to just
# before the first base with: score < q, and prints
# reads whose length is no less than: length > l.

from optparse import OptionParser
import sys, os.path

def checkargs(opts, inFastq, outFastq):
    if len(inFastq) != len(outFastq):
        raise Exception("input the same number of input/output files")
    if len(inFastq) > 1 and len(inFastq) != 2:
        raise Exception("input two intput/output files for paired-end reads")
    for out in outFastq:
        if os.path.exists(out):
            raise Exception("output " + out + " exists")
    if opts.fqtype.lower() == "sanger":
        return 33
    if opts.fqtype.lower() == "illumina":
        return 64
    raise Exception("input correct fastq type")

def dotrim(entry, cutoff, base, minlen):
    numN = 0
    for i in range(len(entry[3])):
        qscore = ord(entry[3][i]) - base
        if cutoff >= qscore:
            seq  = entry[1][:i]
            qual = entry[3][:i]
            if len(seq) == numN:
                return []
            if len(seq) < minlen:
                return []
            part = "0,"+str(i)
            return [entry[0]+":"+part, seq, entry[2], qual]
        if entry[1][i].upper() == 'N':
            numN += 1
    seq = entry[1]
    qual = entry[3]
    part = "0,"+str(len(seq))
    return [entry[0]+":"+part, seq, entry[2], qual]

def readFastq(inFile):
    entry = []
    lcount = 0
    for line in open(inFile):
        line = line.rstrip("\n")
        if line == "": continue;
        lcount += 1
        entry.append(line)
        if lcount == 4:
            yield entry
            lcount = 0
            entry = []

def simpleTrim(opts, args):
    inFastq = args[1].split(",")
    outFastq = args[0].split(",")
    base = checkargs(opts, inFastq, outFastq)
    if len(inFastq) == 2:
        out1 = open(outFastq[0], "w")
        out2 = open(outFastq[1], "w")
        if opts.single_pair == True:
            sg1 = open(outFastq[0]+".single", "w")
            sg2 = open(outFastq[1]+".single", "w")
        for r1, r2 in zip(readFastq(inFastq[0]), readFastq(inFastq[1])):
            r1 = dotrim(r1, opts.cutoff, base, opts.minlen)
            r2 = dotrim(r2, opts.cutoff, base, opts.minlen)
            if r1 != [] and r2 != []:
                out1.write("\n".join(r1)+"\n")
                out2.write("\n".join(r2)+"\n")
            else:
                if opts.single_pair == True:
                    if r1 == [] and r2 != []: sg2.write("\n".join(r2)+"\n");
                    if r2 == [] and r1 != []: sg1.write("\n".join(r1)+"\n");
        out1.close()
        out2.close()
        if opts.single_pair == True:
            sg1.close()
            sg2.close()
    else:
        out = open(outFastq[0], "w")
        for r in readFastq(inFastq[0]):
            r = dotrim(r, opts.cutoff, base, opts.minlen)
            if r != []:
                out.write("\n".join(r)+"\n")
        out.close()

if __name__ == "__main__":
    prog = os.path.basename(sys.argv[0])
    usage = "%prog [options] trimmedOut(s) readFastq(s)"
    description = '''Trim the 3'end of reads to just before the first base below threshould.
For paired-end reads, input two comma-separated fastq files. Trimmed paired-end reads
will be "synchronized" (i.e. all trimmed reads have their pairs).'''

    op = OptionParser(usage=usage, description=description)

    op.add_option("-q","--quality-cutoff",
                  dest="cutoff",
                  type="int",
                  action="store",
                  default=3,
                  help="PHRED quality score cutoff",
                  metavar="VALUE")
    op.add_option("-m", "--min-length",
                  dest="minlen",
                  type="int",
                  action="store",
                  default=30,
                  help="Minimum length of trimmed reads (default=%default)",
                  metavar="VALUE")
    op.add_option("-t", "--fastq-type",
                  dest="fqtype",
                  type="string",
                  action="store",
                  default="sanger",
                  help="PHRED quality score format: 'sanger', 'illumina' (default='%default')",
                  metavar="TYPE")
    op.add_option("-a", "--single-pair",
                  dest="single_pair",
                  action="store_true",
                  default=False,
                  help="Print single paired-end reads (i.e. reads which lose their pairs)")

    (opts, args) = op.parse_args()

    if len(args) != 2:
        raise Exception("specify input/output file names")

    try:
        simpleTrim(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        sys.exit(prog + ": error: " + str(e))

