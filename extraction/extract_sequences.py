#!/usr/bin/env python3

def strip_ambiguities(alignment):
    for sidx in range(len(alignment)):
        seq = alignment[sidx].seq.tomutable()
        for i in range(len(seq)):
            if seq[i] not in "GCAT":
                seq[i] = '-'
        alignment[sidx].seq = seq.toseq()

    return alignment

# Extract sequences

from Bio import Align, AlignIO
from datetime import datetime, timedelta

mixed_alignment = AlignIO.read("Makona_1610_genomes_2016-06-23.fasta", "fasta")

def formatDate(d):
    return d.year + (d - datetime(d.year,1,1)).days/365.0

def getDateStr(seq):
    return seq.id.split("|")[5]

def getLocationStr(seq):
    return seq.id.split("|")[4]

locations = set()

for sidx in range(len(mixed_alignment)):
    seq = mixed_alignment[sidx]

    dateString = getDateStr(seq)
    dt = datetime.strptime(dateString, "%Y-%m-%d")

    seq.annotations["date_string"] = dateString
    seq.annotations["date_numeric"] = formatDate(dt)
    seq.annotations["week"] = dt.isocalendar()[1]

    location = getLocationStr(seq)
    locations.add(location)

for location in locations:

    alignment = Align.MultipleSeqAlignment(
        filter(lambda x: getLocationStr(x) == location, mixed_alignment))

    print("Extracting " + location + " (" + str(len(alignment)) + ") ...")

    strip_ambiguities(alignment)

    # Write fasta output
    AlignIO.write(alignment, "EBOV_" + location + ".fasta", "fasta")
