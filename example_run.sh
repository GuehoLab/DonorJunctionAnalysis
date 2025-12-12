#!/bin/bash


donor_junction_analysis \
    example.fastq.gz \
    output \
    example \
    CCCAACCCCGTGGATGCATTAAGCTGGTCATTGCGGTCTCATTGGTGTACGGTA \
    /mnt/d/Projects/Genome/hg38/hg38_clean \
    chr19:55115754 -t 20
