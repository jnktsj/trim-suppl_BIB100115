#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Junko Tsuji

# This script performs dynamic trimming, the method
# which trims both ends with low quality scores from
# reads.

from optparse import OptionParser
import sys, os.path, re

nPat = re.compile('N', re.IGNORECASE)

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
    rawseq = entry[1]
    rawqual = entry[3]
    lastPos = 0
    beg, end = 0, -1
    for i in range(len(rawqual)):
        qscore = ord(rawqual[i]) - base
        if cutoff >= qscore:
            if len(rawseq[lastPos:i]) > (end-beg):
                beg = lastPos
                end = i
            lastPos = i + 1
    if end == -1:
        end = len(rawseq)
    seq = rawseq[beg:end]
    qual = rawqual[beg:end]
    part = str(beg)+","+str(end)
    L = len(seq)
    if len(nPat.findall(seq)) != L and L >= minlen:
        return [entry[0]+":"+part, seq, entry[2], qual]
    else:
        return []

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

def dynamicTrim(opts, args):
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
    description = '''Trim both ends of reads and take the longest stretch of bases.
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
        dynamicTrim(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        sys.exit(prog + ": error: " + str(e))

