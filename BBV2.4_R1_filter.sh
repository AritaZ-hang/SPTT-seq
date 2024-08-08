#BSUB -q normal
#BSUB -J BC_filter
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=10]"
#BSUB -n 10

cd `pwd`
tmpdir=tmp
mkdir -p $tmpdir
sample_name=$(basename `pwd`)
bbmap_root=/share/home/guoguoji/tools/bbmap/

workdir=`pwd`
echo $workdir

linker1_seq=CGACTCACTACACGT
linker2_seq=TCGCTGACACGATCG

### filter R1 linker ###
${bbmap_root}/bbduk2.sh in=H_R1.fastq.gz in2=H_R2.fastq.gz outm=H_R1.linker1.fastq outm2=H_R2.linker1.fastq fliteral=$linker1_seq k=15 hdist=2

${bbmap_root}/bbduk2.sh in=H_R1.linker1.fastq in2=H_R2.linker1.fastq outm=H_R1.linker2.fastq outm2=H_R2.linker2.fastq fliteral=$linker2_seq k=15 hdist=2 && rm H_R1.linker1.fastq H_R2.linker1.fastq

gzip H_R1.linker2.fastq
gzip H_R2.linker2.fastq

# filter polyT
cutadapt \
--minimum-length 30 \
-a T{10} \
--pair-filter=any \
-o H_R1.use.fastq.gz \
-p H_R2.use.fastq.gz \
--cores 12 H_R1.linker2.fastq.gz H_R2.linker2.fastq.gz > polyT_trim_report.txt
