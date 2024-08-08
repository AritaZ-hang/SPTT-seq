#!/bin/bash


workdir=/media/ggj/Guo-4T-G/fastq_encoding/useData/0804/0618-02/
script=/media/ggj/Guo-4T-G/ST_benchmark/Upstream_codes/SPTT_seq/BBV_pipe/BBV2.4/New_Struct/BBV2.4_mergeFastq.py
paired=/media/ggj/Guo-4T-G/ST_benchmark/Upstream_codes/SPTT_seq/BBV_pipe/BBV2.4/New_Struct/fastqCombinedPairedEnd.py

cd $workdir

python -u $script -bc1 H_R1.bc1.fastq.gz -bc2 H_R1.bc2.fastq.gz -bc3 H_R1.bc3UMI.fastq.gz -technique bbv3.1 -r2 H_R2.use.fastq.gz -out_fq_1 H_R1.final.fastq.gz -out_fq_2 H_R2.final.fastq.gz

python -u $paired H_R1.final.fastq.gz H_R2.final.fastq.gz