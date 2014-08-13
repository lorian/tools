#!/bin/bash
# convert all SAMs to BAMs
for file in *.SAM
    do
    echo $file
    samtools view -hSb $file > ${file/.SAM/.bam}
    done
