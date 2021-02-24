#!/usr/bin/bash
#PBS -j oe
#
#PBS -o /cluster/home/lihe/My_project/BSSeq/log
#PBS -e /cluster/home/lihe/My_project/BSSeq/log
#PBS -l nodes=1:ppn=8

SA=$1
source /cluster/home/lihe/.bashrc
#1,$ s/\${ARRAY\[\$SGE_TASK_ID\]\}/${SA}/g


DIR=/cluster/home/lihe/My_project/BSSeq
RE=$DIR/results
DATA=$DIR/data
TMPD=$DIR/tmp_data
DOC=$DIR/doc
SRC=$DIR/src

REF=/cluster/group/zhujiankang/lihe/Genome/ATH/DNA/arab10.C.fa
REF1=/cluster/group/zhujiankang/lihe/Genome/ATH/DNA/arab10.C.mod.fa

mkdir -p $TMPD/${SA}
/cluster/home/lihe/tools/bsmap-2.90/bsmap -a $TMPD/clean_data/${SA}_Rclean_R1.fastq -b $TMPD/clean_data/${SA}_Rclean_R2.fastq -d $REF -o $TMPD/${SA}/${SA}_bsmap.bsp  -v 2 -p 8 -S 1 2>$TMPD/${SA}/.${SA}_1  ### allow 2 mismatch random report one of the mulitple mappping

/cluster/home/lihe//tools/bsmap-2.90/methratio.py  -r -z -p -d $REF -m 1 $TMPD/${SA}/${SA}_bsmap.bsp -o $TMPD/${SA}/${SA}_bsmap_report.txt  2> $TMPD/${SA}/.${SA}_2 #### remove duplicated , only use the proper mapping pairs
perl $SRC/methyratio_report_to_bed_file.pl $TMPD/${SA}/${SA}_bsmap_report.txt > $TMPD/${SA}/${SA}_methy.bed.original

perl  $SRC/filter_original_bed.pl $TMPD/${SA}/${SA}_methy.bed.original > $TMPD/${SA}/${SA}_methy.bed   ### depth filtering
cd $TMPD/${SA}

# 50 bp bins 
intersectBed -a  ${SA}_methy.bed -b $DOC/TAIR10.bin50bp.bed  -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8,4 -o sum,sum,count|awk '{OFS="\t";print $2,$1,$3,$4,$5}' |perl $SRC/add_class_groupBy.weight.pl|awk -v a=${SA} '{print $0"\t"a}' >  ${SA}_bin50bp.bed                                        

#Gene_body
intersectBed -a  ${SA}_methy.bed -b $DOC/TAIR10.C.gene.bed -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.C.gene.bed >  ${SA}_gene_ML.out
# Protein coding gene 
perl   ~/PC/code/join_A_B_byteB.pl $DOC/Protein_coding_gene.ID ${SA}_gene_ML.out > ${SA}_protein_coding_gene_ML.out

#Gene_up_2K
intersectBed -a  ${SA}_methy.bed -b $DOC/TAIR10.C.gene_up2k.bed -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.C.gene_up2k.bed  >  ${SA}_gene_up2k_ML.out

###Gene_down_2K
intersectBed -a  ${SA}_methy.bed -b $DOC/TAIR10.C.gene_down2k.bed -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.C.gene_down2k.bed  >  ${SA}_gene_down2k_ML.out

## UP and Down 2k
intersectBed -a  ${SA}_methy.bed -b /cluster/group/zhujiankang/lihe/Genome/ATH//annotation/TAIR10.C.gene.nearBy.bed -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl $SRC/add_class_groupBy.pl /cluster/group/zhujiankang/lihe/Genome/ATH//annotation/TAIR10.C.gene.nearBy.bed  >  ${SA}_gene_nearBy_ML.out

#/oldgroupshare/bioinformatics/zhaocheng/Genome/ATH/TE/TE_bedlike.bed
intersectBed -a  ${SA}_methy.bed -b <(cut -f 1-6 $DOC/TE_bedlike.bed) -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl  $SRC/add_class_groupBy.pl $DOC/TE_bedlike.bed >  ${SA}_TE.out

# IG region
intersectBed -a  ${SA}_methy.bed -b $DOC/TAIR10.IG.bed -wa -wb |sort -k5,5 -k12,12 |groupBy  -g 5,12 -c 7,8 -o sum,sum|awk '{OFS="\t";print $2,$1,$3,$4}' |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.IG.bed >  ${SA}_IG.out

# Gene 60 bp bin
intersectBed -a <(awk -F "\t" '{OFS="\t";print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TAIR10.C.gene.bed|perl $SRC/UD_modify.pl) -b  <(awk -F "\t" '{OFS="\t";print $1,$2,$3,$4"|"$5"|"$7"|"$8,$5,$6}' ${SA}_methy.bed) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g'|awk -F "\t" '{print $1"_"$6"_"$7"\t"$3"\t"$4"\t"$5}' |sort -k1,1 -k2,2 |groupBy -g 1,2 -c 3,4 -o sum,sum |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.C.gene.nearBy.60bin.out >  ${SA}_gene.nearBy.60bin.out


### TE 60 bp bin
intersectBed -a <(awk -F "\t" '{OFS="\t";if (($3-$2) >99) print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TE_bedlike.bed|perl $SRC/UD_modify.pl) -b  <(awk -F "\t" '{OFS="\t";print $1,$2,$3,$4"|"$5"|"$7"|"$8,$5,$6}' ${SA}_methy.bed) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g'|awk -F "\t" '{print $1"_"$6"_"$7"\t"$3"\t"$4"\t"$5}' |sort -k1,1 -k2,2 |groupBy -g 1,2 -c 3,4 -o sum,sum |perl $SRC/add_class_groupBy.pl $DOC/TAIR10.TE.nearBy.60bin.out >  ${SA}_TE.nearBy.60bin.out

#
rm ${SA}_bsmap.bsp ${SA}_bsmap_report.txt



for ITEM in  ${SA}_TE.out ${SA}_protein_coding_gene_ML.out ${SA}_IG.out
do
	echo "${SA}	`sort -k2,2  ${ITEM} |groupBy  -g 2 -c 3,4 -o sum,sum |tr '\n' '\t'|awk '{OFS="\t";print $2/($2+$3),$5/($5+$6),$8/($8+$9),$11/($11+$12)}'`" > ${ITEM}.me.part
done

# end of job script

