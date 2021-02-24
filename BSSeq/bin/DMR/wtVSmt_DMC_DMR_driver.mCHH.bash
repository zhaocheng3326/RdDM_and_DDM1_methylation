#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq
TMPD=$DIR/tmp_data
SRC=$DIR/src/
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin


SA=$1
SB=$2
cd $TMPD


mkdir -p DMR/${SB}_VS_${SA}_CHH

cut -f 1-3,5,6 ${SA}/${SA}_methy.bed ${SB}/${SB}_methy.bed|grep "CHH" |sort |uniq -d  > DMR/${SB}_VS_${SA}_CHH/wt_mt.step1.bed

perl $SRC/Step1_in_two_library.pl ${SA}/${SA}_methy.bed DMR/${SB}_VS_${SA}_CHH/wt_mt.step1.bed  > DMR/${SB}_VS_${SA}_CHH/wt_methy.bed.step1
perl $SRC/Step1_in_two_library.pl ${SB}/${SB}_methy.bed DMR/${SB}_VS_${SA}_CHH/wt_mt.step1.bed  > DMR/${SB}_VS_${SA}_CHH/mt_methy.bed.step1

cd  $TMPD/DMR/${SB}_VS_${SA}_CHH
paste wt_methy.bed.step1 mt_methy.bed.step1|cut -f 1-3,5-8,15-16 > wt_mt.DMC_input.step2
perl $SRC/Step2_DMC_cal.pl wt_mt.DMC_input.step2 > wt_mt.DMC_output.step3
awk '{if (($11 == "Hyper")&& ($10 <0.01)) print $0}' wt_mt.DMC_output.step3 > wt_mt.DMC_output.step3.DMC.Hyper
awk '{if (($11 == "Hypo")&& ($10 <0.01)) print $0}' wt_mt.DMC_output.step3 > wt_mt.DMC_output.step3.DMC.Hypo
cat wt_mt.DMC_output.step3.DMC.Hyper wt_mt.DMC_output.step3.DMC.Hypo |sortBed > wt_mt.DMC_output.step3.DMC


intersectBed -a <(intersectBed -a $DOC/TAIR10.windows.mCHH.sel.bed -b wt_mt.DMC_output.step3.DMC -wa|sort|uniq ) -b wt_methy.bed.step1 -wa -wb|cut -f 1-4,9,11,12 > wt_methy.bed.DMR.step4
intersectBed -a <(intersectBed -a $DOC/TAIR10.windows.mCHH.sel.bed -b wt_mt.DMC_output.step3.DMC -wa|sort|uniq ) -b mt_methy.bed.step1 -wa -wb|cut -f 1-4,9,11,12 > mt_methy.bed.DMR.step4

sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5  wt_methy.bed.DMR.step4 |groupBy -g 1,2,3,4,5 -c  6,7 -o sum,sum > wt_methy.bed.DMR.step5  ### split
sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5  mt_methy.bed.DMR.step4 |groupBy -g 1,2,3,4,5 -c  6,7 -o sum,sum > mt_methy.bed.DMR.step5  ### split
#sort -k 1,1 -k2,2n -k3,3n -k4,4 wt_methy.bed.DMR.step5 |groupBy -g 1,2,3,4 -c  6,7 -o sum,sum |awk '{OFS="\t";print $1,$2,$3,$4,"CX",$5,$6}' > wt_methy.bed.DMR.CX.step5
#sort -k 1,1 -k2,2n -k3,3n -k4,4 mt_methy.bed.DMR.step5 |groupBy -g 1,2,3,4 -c  6,7 -o sum,sum |awk '{OFS="\t";print $1,$2,$3,$4,"CX",$5,$6}' > mt_methy.bed.DMR.CX.step5

#cat wt_methy.bed.DMR.step5 wt_methy.bed.DMR.CX.step5 |sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5 > wt_methy.bed.DMR.step5B
#cat mt_methy.bed.DMR.step5 mt_methy.bed.DMR.CX.step5 |sort -k 1,1 -k 2,2n -k 3,3n -k 4,4 -k 5,5 > mt_methy.bed.DMR.step5B

paste wt_methy.bed.DMR.step5 mt_methy.bed.DMR.step5|cut -f 1-7,13,14 |grep "CHH" > wt_mt_methy.bed.DMR_CHH_input.step6





perl $SRC/Step7_DMR_cal_abs.pl wt_mt_methy.bed.DMR_CHH_input.step6 > Step7.temp                                                                                                                           
/cluster/tools/R-3.3.3/bin/Rscript $SRC/Step7_DMR_cal.abs.R Step7.temp CHH.abs.step7

perl $SRC/Step7_DMR_cal.pl wt_mt_methy.bed.DMR_CHH_input.step6 > Step7.temp 

/cluster/tools/R-3.3.3/bin/Rscript $SRC/Step7_DMR_cal.pl.R Step7.temp CHH.logFC.step7

intersectBed -a <(awk '{if (($12 < 1.1) &&  ($10 > 0)) print $0}' CHH.abs.step7) -b wt_mt.DMC_output.step3.DMC.Hyper -wa -wb |cut -f 1-12 |perl $SRC/count_repeat.pl  > CHH.abs.step8
intersectBed -a <(awk '{if (($12 < 1.1) &&  ($10 < 0)) print $0}' CHH.abs.step7) -b wt_mt.DMC_output.step3.DMC.Hypo -wa -wb |cut -f 1-12 |perl $SRC/count_repeat.pl  >> CHH.abs.step8
intersectBed -a <(awk '{if (($12 < 1.1) &&  ($10 > 0)) print $0}' CHH.logFC.step7) -b wt_mt.DMC_output.step3.DMC.Hyper -wa -wb |cut -f 1-12 |perl $SRC/count_repeat.pl  > CHH.logFC.step8
intersectBed -a <(awk '{if (($12 < 1.1) &&  ($10 < 0)) print $0}' CHH.logFC.step7) -b wt_mt.DMC_output.step3.DMC.Hypo -wa -wb |cut -f 1-12 |perl $SRC/count_repeat.pl  >> CHH.logFC.step8	
	
	

