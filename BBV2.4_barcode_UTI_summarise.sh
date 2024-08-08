#BSUB -J UTI_summarises
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 36

summariseScript=/share/home/guoguoji/RAWDATA/index_microscopy/scripts/BBV2.4_barcode_UTI_summarise.py

workdir=/share/home/guoguoji/RAWDATA/index_microscopy/0717/0705-01/WithLinker/

in_fq_1=H_R1_filtered.fastq.gz
in_fq_2=H_R2_filtered.fastq.gz

python -u $summariseScript -in_fq_1 $in_fq_1 -in_fq_2 $in_fq_2 -workdir $workdir