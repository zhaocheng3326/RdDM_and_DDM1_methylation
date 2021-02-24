#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq/
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin




P1=( head abs)  ##### cutoff
P2=( head 3 4 5 )    ##### DMC
P3=( head 0.1 0.08 0.05)

P1=( head logFC )
P2=( head 3 4 5 )
P3=( head 0.58 1 )


cd $DIR


ARRAY=( head Col Col Col  Col  Col  d4 ddm1 ddm1 pol4 pol4 ddm1 pol4 pol5 pol5 ddm1 Col  Col      ddm1     ddm1 cmt2     d4       d5 )
BRRAY=( head d4  d5  ddm1 pol4 pol5 d5 d4   d5   d4   d5   pol4 pol5 d4   d5   pol5 cmt2 ddm1cmt2 ddm1cmt2 cmt2 ddm1cmt2 ddm1cmt2 ddm1cmt2 )

DF=logFC
for DMC in  4
do
	for RA in 1.58
	do
		for I in {1..15}
		do
			mkdir -p $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}
			cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}

			python3 $SRC/select_windows.py $TMPD/DMR/${BRRAY[$I]}_merge_VS_${ARRAY[$I]}_merge_CHH/CHH.${DF}.step8 ${DMC} ${RA} 0.05  > rep_conserve.out #

		done
	done
done


for DMC in  4
do
	for RA in 1.58
	do
		for I in 16 17 18 19 21 22
		do
			mkdir -p $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}
			cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}
			python3 $SRC/select_windows.py $TMPD/DMR/${BRRAY[$I]}_VS_${ARRAY[$I]}_merge_CHH/CHH.${DF}.step8 ${DMC} ${RA} 0.05  > rep_conserve.out
		done
		cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/
	done
done

for DMC in  4
do
	for RA in 1.58
	do
		for I in 20
		do
			mkdir -p $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}
			cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/${BRRAY[$I]}_VS_${ARRAY[$I]}

			python3 $SRC/select_windows.py $TMPD/DMR/${BRRAY[$I]}_VS_${ARRAY[$I]}_CHH/CHH.${DF}.step8 ${DMC} ${RA} 0.05  > rep_conserve.out
		done
		cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/

	done
done

cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/

mkdir -p temp
cat pol4_VS_Col/rep_conserve.out pol5_VS_Col/rep_conserve.out |grep "Hypo" |perl ~/PC/code/count_repeat.pl|awk '{if ($3 >=2) print $1"\t"$2}' >  temp/Pol_VS_Col.Hypo.out #### classic

cat pol4_VS_Col/rep_conserve.out pol5_VS_Col/rep_conserve.out ddm1_VS_Col/rep_conserve.out  |grep -v "Hypo" |cut -f 1  |perl ~/PC/code/count_repeat.pl|awk '{if ($2 ==3) print $1"\tNot_change"}' > temp/single_mutant.NC.out

cat d4_VS_Col/rep_conserve.out d5_VS_Col/rep_conserve.out|perl ~/PC/code/count_repeat.pl|grep "Hypo"|awk '{if ($3>=2) print $1"\t"$2}' >  temp/Double.Hypo.out

cat temp/Double.Hypo.out <(cat temp/single_mutant.NC.out | sed -e 's/Not_change/Hypo/' ) |perl ~/PC/code/count_repeat.pl|awk '{if ($3==2) print $1"\t"$2}' > temp/Only.Double.Hypo.out

cat d4_VS_ddm1/rep_conserve.out d5_VS_ddm1/rep_conserve.out|grep  "Hypo"| perl ~/PC/code/count_repeat.pl|awk '{if ($3 ==2) print $1}' > temp/d45_vs_ddm1.Hypo.out #tag1
cat ddm1cmt2_VS_d4/rep_conserve.out ddm1cmt2_VS_d5/rep_conserve.out|grep  "Hyper"| perl ~/PC/code/count_repeat.pl|awk '{if ($3 ==2) print $1}' > temp/d45_vs_ddm1cmt2.Hypo.out #tag2

cat cmt2_VS_Col/rep_conserve.out |grep "Hypo"  > temp/cmt2.Hypo.out
cat ddm1cmt2_VS_Col/rep_conserve.out |grep "Hypo"  > temp/ddm1cmt2.Hypo.out

bash $BIN/suit/over_setdiff.bash setdiff1 temp/Pol_VS_Col.Hypo.out  temp/cmt2.Hypo.out |cut -f 4 > temp/Pol-cmt2.Hypo.out
bash $BIN/suit/over_setdiff.bash setdiff1 temp/cmt2.Hypo.out  temp/Pol_VS_Col.Hypo.out   |cut -f 4 > temp/cmt2-Pol.Hypo.out
bash $BIN/suit/over_setdiff.bash overlap temp/cmt2.Hypo.out  temp/Pol_VS_Col.Hypo.out  > temp/cmt2_ov_Pol.Hypo.out
bash $BIN/suit/over_setdiff.bash overlap temp/Double.Hypo.out temp/cmt2-Pol.Hypo.out  > temp/d45_ov___p__cmt2-Pol__p__.Hypo.out

#way1
bash $BIN/suit/over_setdiff.bash overlap temp/d45_ov___p__cmt2-Pol__p__.Hypo.out temp/d45_vs_ddm1.Hypo.out > temp/d45_ov___p__cmt2-Pol__p__.Hypo._ov_tag1.out
bash $BIN/suit/over_setdiff.bash overlap temp/d45_ov___p__cmt2-Pol__p__.Hypo._ov_tag1.out temp/d45_vs_ddm1cmt2.Hypo.out> temp/d45_ov___p__cmt2-Pol__p__.Hypo._ov_tag1_ov_tag2.out #tag3, cand C3
bash $BIN/suit/over_setdiff.bash  setdiff1 temp/cmt2-Pol.Hypo.out temp/d45_ov___p__cmt2-Pol__p__.Hypo._ov_tag1_ov_tag2.out |cut -f 4 > temp/cmt2-Pol.Hypo.-tag3.out # cand C2
#way2
bash $BIN/suit/over_setdiff.bash  setdiff1 temp/d45_ov___p__cmt2-Pol__p__.Hypo.out  temp/ddm1cmt2.Hypo.out|cut -f 4 > temp/d45_ov___p__cmt2-Pol__p__-ddm1cmt2.Hypo.out #tag4 #cand C3
bash $BIN/suit/over_setdiff.bash  setdiff1 temp/cmt2-Pol.Hypo.out temp/d45_ov___p__cmt2-Pol__p__-ddm1cmt2.Hypo.out |cut -f 4 > temp/cmt2-Pol.Hypo.-tag4.out #cand C2


#choose way2
cut -f 1 temp/Pol-cmt2.Hypo.out > temp/C1.out
cut -f 1 temp/cmt2-Pol.Hypo.-tag4.out > temp/C2.out
cut -f 1 temp/d45_ov___p__cmt2-Pol__p__-ddm1cmt2.Hypo.out > temp/C3.out
cut -f 1 temp/cmt2_ov_Pol.Hypo.out > temp/C4.out

bash $BIN/suit/over_setdiff.bash report temp/C1.out > temp/C1.temp.bed
bash $BIN/suit/over_setdiff.bash report temp/C2.out > temp/C2.temp.bed
bash $BIN/suit/over_setdiff.bash report temp/C3.out > temp/C3.temp.bed
bash $BIN/suit/over_setdiff.bash report temp/C4.out > temp/C4.temp.bed

cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/
echo -n "" > total.analysis.200version.bed
for ITEM in C1 C2 C3 C4
do
	sortBed -i temp/${ITEM}.temp.bed |mergeBed|awk -v a=${ITEM} '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3"\t"a}' >> total.analysis.200version.bed
done


cd $TMPD/DMR_CHH_${DF}_${DMC}_${RA}/
intersectBed -a <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 $TMPD/Col_meCHH.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk '{OFS="\t";print $4,$2,$3}') -b total.analysis.200version.bed -f 0.99 -wa -wb |sort -k7,7 -k8,8 -k1,1 |groupBy -g 1,7,8 -c 2,3 -o min,max |awk '{OFS="\t";print $1,$4,$5,$1"|"$4"|"$5"|"$3}'  > total.analysis.QS.temp

echo -n "" > total.analysis.QS.bed
for CLA in C1 C2 C3 C4
do
	grep ${CLA} total.analysis.QS.temp |sortBed |mergeBed |awk -v a=${CLA}  '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3"|"a}' |awk '{if (($3-$2)>= 50) print $0}' >> total.analysis.QS.bed
done

cd /cluster/home/lihe/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58/
grep "C1" total.analysis.QS.bed |sed -e 's/|C1/_C1/'> C1.bed
grep "C2" total.analysis.QS.bed |sed -e 's/|C2/_C2/'> C2.bed
grep "C3" total.analysis.QS.bed |sed -e 's/|C3/_C3/'> C3.bed
grep "C4" total.analysis.QS.bed |sed -e 's/|C4/_C4/'> C4.bed
