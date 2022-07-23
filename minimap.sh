#!/bin/bash
echo "export CURRTIME=$(date +%s)" > /tmp/curr_time
time ./minimap2-2.24_x64-linux/minimap2 -ax map-ont --MD -t 14 chr17.fa n_guppy.fq | samtools view -F 0x904 -hb -@6 | samtools sort -@6 | samtools view -hb -@6 > all.bam
time samtools index -@14 all.bam
rm samtools.*
