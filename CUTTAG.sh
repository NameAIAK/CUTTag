cd /data/CUTTag/n87
mkdir 1_rawdata 2_clean 3_sam 3_bam 4_rmdup 5_bw 6_callpeak

#  复制数据到1_rawdata；数据命名：N87_1_R1.fq.gz    N87_1_R2.fq.gz  sample_R1.fq.gz sample_R2.fq.gz
# 1. clean data : 
conda activate RNA-seq
ls 1_rawdata/*_R1.fq.gz |while read id ; do trim_galore -j 50 -q 30 --phred33 --length 120 -e 0.1  --paired -o 2_clean/ $id ${id%_R1.*gz}_R2.fq.gz;done

# 2.map: 
#小鼠数据库/data/RNA/reference/bowtie2_mouse/mm10；人源数据库/data/RNA/reference/bowtie2_human/GRCh38/GRCh38
cd 2_clean
#人人人人人人人人人人人人人人人人人人人人
ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /data/RNA/reference/bowtie2_human/GRCh38/GRCh38  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done
#小鼠小鼠小鼠小鼠小鼠小鼠小鼠小鼠小鼠小鼠
ls *R1*.fq.gz  | while read id ; do bowtie2  -p 50  -x /data/RNA/reference/bowtie2_mouse/mm10  --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -1 ${id%_R1*}_R1_val_1.fq.gz -2 ${id%_R1*}_R2_val_2.fq.gz  -S ../3_sam/${id%_R1*}.sam;done

# 3.samtobam:
cd ../3_sam
ls *.sam |while read id ;do samtools view -@ 50 -bF 12 -q 30 -S $id > ../3_bam/${id%.sam}.map.q30.F12.bam;done

# 3-2.remove duplication: 
cd ../3_bam
ls *.bam |while read id ;do /home/star/anaconda3/envs/sambamba/bin/sambamba markdup -r $id ../4_rmdup/${id%.bam}.rmdup.bam ;done

# 4 sort bam and make index
cd ../4_rmdup
ls *.rmdup.bam |while read id ;do samtools sort -@ 50 $id -o ${id%.rmdup.bam}.rmdup.sorted.bam;done
ls *.rmdup.sorted.bam |while read id ;do samtools index $id;done

# 5.bam to bw: 
ls *.rmdup.sorted.bam |while read id ;do  bamCoverage -b $id --normalizeUsing RPKM -p 50 -o ../5_bw/${id%.map.q30.F12.rmdup.sorted.bam}.trimgalore.map.q30.F12.rmdup.rpkm.bw ; done

# 6.call peak(adjust for transcrip factor/ Histone modification type)：
conda activate macs2
sample1='N87_1'
sample2='N87_2'
pre='N87'
macs2  callpeak -t $sample1'.map.q30.F12.rmdup.sorted.bam' $sample2'.map.q30.F12.rmdup.sorted.bam' --bdg -p 1e-5 -g hs -n '../6_callpeak/'$pre'.peak'

# 7.visualization
conda activate ChIP-seq
cd ../6_callpeak
computeMatrix  reference-point -p 20 -R $pre'.peak_summits.bed' -a 3000 -b 3000 -S  '../5_bw/'$sample1'.trimgalore.map.q30.F12.rmdup.rpkm.bw' '../5_bw/'$sample2'.trimgalore.map.q30.F12.rmdup.rpkm.bw' --skipZeros  -out ./$pre'.computeMatrix.gz'
#通过此图可以查看样本测序情况
plotHeatmap -m $pre'.computeMatrix.gz' -o $pre'.computeMatrix.gz.pdf' --colorMap RdBu --zMin -3 --zMax 3

# 8.annotation
conda activate macs2
annotatePeaks.pl $pre'.peak_summits.bed' hg38 > $pre'_peak.anno.txt'
tail +2 '../6_callpeak/'$pre'_peak.anno.txt' |cut -f 2,3,4 > '../6_callpeak/'$pre'_peak.txt'

# 9.bam2tdf
conda activate RNA-seq
igvtools count ../4_rmdup/$sample1'.map.q30.F12.rmdup.sorted.bam '$pre'_1.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1
igvtools count ../4_rmdup/$sample2'.map.q30.F12.rmdup.sorted.bam '$pre'_2.tdf' hg38.chrom.sizes 1>>bam2tdf.log 2>&1

