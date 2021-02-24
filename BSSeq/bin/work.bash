#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq/
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin
LOG=$DIR/log

cd $DIR
### make  the meta.table 
echo -e "Sample\tGSample\tReplicate" > $DOC/BSSeq.metadata.txt
for ITEM  in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 ddm1cmt2 ddm1cmt3 ddm1drm1 drm12 nrpd1cmt2 nrpd1cmt3 pol4_1 pol4_2 pol5_1 pol5_2
do
	echo ${ITEM}|grep "_"|awk '{print $0"\t"$0}'|sed -e 's/_/\t/'|awk -F "\t" '{OFS="\t";print $3,$1,"BRep_"$2}' >> $DOC/BSSeq.metadata.txt
	echo ${ITEM}|grep -v "_"|awk -F "\t" '{OFS="\t";print $1,$1,"BRep_0"}' >> $DOC/BSSeq.metadata.txt
done
#### QC
#for ITEM in cmt2_R Col_1_R Col_2_R d4_1_R d4_2_R d5_1_R d5_2_R ddm1_1_R ddm1_2_R ddm1cmt2_R ddm1cmt3_R ddm1drm1_R drm12_R nrpd1cmt2_R nrpd1cmt3_R pol4_1_R pol4_2_R pol5_1_R pol5_2_R
#do
#	sbatch -J ${ITEM}.QC -e ${LOG}/${ITEM}.QC.log -o ${LOG}/${ITEM}.QC.log  --mincpus=2 $BIN/QC/QC.bash ${ITEM} 
#	sleep 10
#done

###BS_map_step1
#for ITEM in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 ddm1cmt2 ddm1cmt3 ddm1drm1 drm12 nrpd1cmt2 nrpd1cmt3 pol4_1 pol4_2 pol5_1 pol5_2
#do
#	sbatch -J ${ITEM}.BS1 -e ${LOG}/${ITEM}.BS1.log -o ${LOG}/${ITEM}.BS1.log  --mincpus=4 $BIN/Mapping/BSMAP_step1.bash ${ITEM} 
#	sleep 10
#done

###BS_map_step2
for ITEM  in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 ddm1cmt2 ddm1cmt3 ddm1drm1 drm12 nrpd1cmt2 nrpd1cmt3 pol4_1 pol4_2 pol5_1 pol5_2
do
	sbatch -J BS2.${ITEM} -e ${LOG}/${ITEM}.BS2.log -o ${LOG}/${ITEM}.BS2.log  --mincpus=4 $BIN/Mapping/BSMAP_step2.bash ${ITEM} 
	sleep 10
done

# merge replicate
ARRAY=( head  Col ddm1 pol4 pol5 d4 d5 )
for n in {1..6}
do
	cd $TMPD
	mkdir -p ${ARRAY[$n]}_merge
	cd $TMPD/${ARRAY[$n]}_1/
	awk -F "\t" '{print>$1}' ${ARRAY[$n]}_1_methy.bed.original
	cd $TMPD/${ARRAY[$n]}_2/
	awk -F "\t" '{print>$1}' ${ARRAY[$n]}_2_methy.bed.original
	cd $TMPD
	echo -n "" > ${ARRAY[$n]}_merge/${ARRAY[$n]}_merge_methy.bed.original 
	for CHR in Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
	do
		cat ${ARRAY[$n]}_1/${CHR} ${ARRAY[$n]}_2/${CHR}|sort -k1,1 -k2,2 -k3,3 -k5,5 -k6,6  |groupBy -g 1,2,3,5,6 -c 7,8 -o sum,sum |awk -F "\t" '{OFS="\t";print $1,$2,$3,$4"."NR,$4,$5,$6,$7}' >> ${ARRAY[$n]}_merge/${ARRAY[$n]}_merge_methy.bed.original
	done
	cd $TMPD
	
	perl  $SRC/filter_original_bed.pl ${ARRAY[$n]}_merge/${ARRAY[$n]}_merge_methy.bed.original > ${ARRAY[$n]}_merge/${ARRAY[$n]}_merge_methy.bed
done
	

ARRAY=( head  Col_merge ddm1_merge pol4_merge pol5_merge cmt2 d4_merge d5_merge ddm1cmt2 )
for m in {1..8}
do
	for n in {1..8}
	do
		if [ $n -gt $m ];then
		echo "${ARRAY[$m]} ${ARRAY[$n]}"
		sbatch -J ${ARRAY[$m]}_${ARRAY[$n]}.DMR -e ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.log -o ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.log --cpus-per-task 2 --mincpus=2 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ${ARRAY[$m]} ${ARRAY[$n]}
		sleep 2		
		fi;
	done
done

# replicate split 
#ARRAY=( head  Col ddm1 pol4 pol5 d4 d5 )
#for m in {1..6}
#do
#	for n in {1..6}
#	do
#		if [ $n -gt $m ];then
#		echo "${ARRAY[$m]} ${ARRAY[$n]}"
#		sbatch -J ${ARRAY[$m]}_${ARRAY[$n]}.rep1.DMR -e ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.rep1.log -o ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ${ARRAY[$m]}_1 ${ARRAY[$n]}_1
#		sleep 5
#		sbatch -J ${ARRAY[$m]}_${ARRAY[$n]}.rep2.DMR -e ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.rep2.log -o ${LOG}/${ARRAY[$m]}_${ARRAY[$n]}.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ${ARRAY[$m]}_2 ${ARRAY[$n]}_2
#		sleep 5
		
#		fi;
#	done
#done

#sbatch -J ddm1cmt2_Col.rep1.DMR -e ${LOG}/ddm1cmt2_Col.rep1.log -o ${LOG}/ddm1cmt2_Col.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash Col_1 ddm1cmt2
#sleep 10
#sbatch -J ddm1cmt2_Col.rep2.DMR -e ${LOG}/ddm1cmt2_Col.rep2.log -o ${LOG}/ddm1cmt2_Col.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash Col_2 ddm1cmt2

#sbatch -J ddm1cmt2_ddm1.rep1.DMR -e ${LOG}/ddm1cmt2_ddm1.rep1.log -o ${LOG}/ddm1cmt2_ddm1.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ddm1_1 ddm1cmt2
#sleep 10
#sbatch -J ddm1cmt2_ddm1.rep2.DMR -e ${LOG}/ddm1cmt2_ddm1.rep2.log -o ${LOG}/ddm1cmt2_ddm1.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ddm1_2 ddm1cmt2

#sbatch -J cmt2_ddm1.rep1.DMR -e ${LOG}/cmt2_ddm1.rep1.log -o ${LOG}/cmt2_ddm1.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ddm1_1 cmt2
#sleep 10
#sbatch -J cmt2_ddm1.rep2.DMR -e ${LOG}/cmt2_ddm1.rep2.log -o ${LOG}/cmt2_ddm1.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ddm1_2 cmt2

#sbatch -J ddm1cmt2_cmt2.DMR -e ${LOG}/ddm1cmt2_cmt2.log -o ${LOG}/ddm1cmt2_cmt2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash ddm1cmt2 cmt2

#sbatch -J ddm1cmt2_d4.rep1.DMR -e ${LOG}/ddm1cmt2_d4.rep1.log -o ${LOG}/ddm1cmt2_d4.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash d4_1 ddm1cmt2
#sbatch -J ddm1cmt2_d4.rep2.DMR -e ${LOG}/ddm1cmt2_d4.rep2.log -o ${LOG}/ddm1cmt2_d4.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash d4_2 ddm1cmt2

#sbatch -J ddm1cmt2_d5.rep1.DMR -e ${LOG}/ddm1cmt2_d5.rep1.log -o ${LOG}/ddm1cmt2_d5.rep1.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash d5_1 ddm1cmt2
#sbatch -J ddm1cmt2_d5.rep2.DMR -e ${LOG}/ddm1cmt2_d5.rep2.log -o ${LOG}/ddm1cmt2_d5.rep2.log --mincpus=4 $BIN/DMR/wtVSmt_DMC_DMR_driver.mCHH.bash d5_2 ddm1cmt2
##







### Expression and correlation
# intersectBed -a TAIR10.bin50bp.bed  -b <(cat TAIR10.C.gene.bed  TE_bedlike.bed |cut -f 1-4 |awk '{OFS="\t";if ($2 <2000) print $1,0,$3,$4;else print $1,$2-2000,$3+2000,$4}' ) -wa -wb |cut -f 4,10 |sort |uniq  > bin50bp.asso.gene
# cut -f 1  ~/My_project/RNASeq/tmp_data/DEG/*DEG|grep -v "Gene"|sort |uniq  >  ~/My_project/RNASeq/tmp_data/DEG/DEG.ID.total
# join -1 2 -2 1 -t$'\t' <(sort -k2,2 $DOC/bin50bp.asso.gene) <(sort -k1,1 ~/My_project/RNASeq/tmp_data/DEG/DEG.ID.total) > $DOC/DEG.sel.bin50bp.asso.gene
echo -e "win\tGene\tType\tC\tT\tWeight\tSample" > DEG.sel.bin50bp.met.out
for ITEM in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 pol4_1 pol4_2
do
	join -1 2 -2 1 -t$'\t' <(sort -k2,2 $DOC/DEG.sel.bin50bp.asso.gene) <(sort -k1,1 $TMPD/${ITEM}/${ITEM}__bin50bp.bed) >> DEG.sel.bin50bp.met.out
done



###
bash $BIN/DMR/parameters_CHH_DMR.analysis.2007.bash
bash $BIN/DMR/DMR_feature.analysis.bash

#bash $BIN/DMR_work_ML.bash
#Rscript ~/My_project/BSSeq/src/Analysis/region_heatmap.R
#Rscript ~/My_project/BSSeq/src/Analysis/region_heatmap1.2.R

