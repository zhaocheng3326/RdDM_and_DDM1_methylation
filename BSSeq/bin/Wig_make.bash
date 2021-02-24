#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq
TMPD=$DIR/tmp_data
SRC=$DIR/src/
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin

REF=/cluster/group/zhujiankang/lihe/Genome/ATH/DNA/arab10.C.fa
CODE3=$SRC/deal_priwig.pl

cd $TMPD
mkdir -p WIG

for ITEM in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 ddm1cmt2 ddm1cmt3 ddm1drm1 drm12 nrpd1cmt2 nrpd1cmt3 pol4_1 pol4_2 pol5_1 pol5_2
do
	cd $TMPD/${ITEM}
	perl $SRC/report_bed_to_priwig.pl ${ITEM}_methy.bed
done



for ITEM in CX CG CHG CHH
do
	for IN in cmt2 Col_1 Col_2 d4_1 d4_2 d5_1 d5_2 ddm1_1 ddm1_2 ddm1cmt2 ddm1cmt3 ddm1drm1 drm12 nrpd1cmt2 nrpd1cmt3 pol4_1 pol4_2 pol5_1 pol5_2
	do
		cd $TMPD/${IN}
		sort -k 1,1 -k 2n,2  ${IN}_${ITEM}.priwig -o ${IN}_${ITEM}.priwig
		perl $CODE3 ${IN}_${ITEM} ${IN}_${ITEM}.priwig|awk '{if ($2 !=0) print $0}'  > ${IN}_${ITEM}.wig
		#/share/apps/prog/IGVTools/igvtools toTDF ${IN}_${ITEM}.wig $TMPD/WIG/${IN}_${ITEM}.wig.tdf  $REF
		rm ${IN}_${ITEM}.priwig 
		mv ${IN}_${ITEM}.wig $TMPD/WIG/
	done
done

for IN in Col d4 d5 ddm1 pol4 pol5
do
  cd $TMPD/${IN}_1
  perl $SRC/report_bed_to_priwig.pl ${IN}_1_methy.bed
  
  cd $TMPD/${IN}_2
  perl $SRC/report_bed_to_priwig.pl ${IN}_2_methy.bed
        
  for ITEM in CX CG CHG CHH
  do
    cd $TMPD/${IN}_1
    sort -k 1,1 -k 2n,2  ${IN}_1_${ITEM}.priwig -o ${IN}_1_${ITEM}.priwig
    cd $TMPD/${IN}_2
    sort -k 1,1 -k 2n,2  ${IN}_2_${ITEM}.priwig -o ${IN}_2_${ITEM}.priwig
  
    cd $TMPD/WIG
    cat $TMPD/${IN}_1/${IN}_1_${ITEM}.priwig $TMPD/${IN}_2/${IN}_2_${ITEM}.priwig |sort -k1,1 -k2,2 |groupBy -g 1,2 -c 3,3 -o sum,count |awk '{OFS="\t";print $1,$2,$3/$4}' > $TMPD/WIG/${IN}_aver_${ITEM}.priwig
    sort -k 1,1 -k 2n,2 $TMPD/WIG/${IN}_aver_${ITEM}.priwig -o $TMPD/WIG/${IN}_aver_${ITEM}.priwig
    
    perl $CODE3 ${IN}_aver_${ITEM} ${IN}_aver_${ITEM}.priwig|awk '{if ($2 !=0) print $0}'  > ${IN}_aver_${ITEM}.wig
  done
done
for ITEM in Col_merge d4_merge d5_merge ddm1_merge pol4_merge pol5_merge
do
  cd $TMPD/${ITEM}
  perl $SRC/report_bed_to_priwig.pl ${ITEM}_methy.bed
done
for ITEM in CX CG CHG CHH
do
    for IN in Col_merge d4_merge d5_merge ddm1_merge pol4_merge pol5_merge
    do
        cd $TMPD/${IN}
        sort -k 1,1 -k 2n,2  ${IN}_${ITEM}.priwig -o ${IN}_${ITEM}.priwig
        perl $CODE3 ${IN}_${ITEM} ${IN}_${ITEM}.priwig|awk '{if ($2 !=0) print $0}'  > ${IN}_${ITEM}.wig
        #/share/apps/prog/IGVTools/igvtools toTDF ${IN}_${ITEM}.wig $TMPD/WIG/${IN}_${ITEM}.wig.tdf  $REF
        rm ${IN}_${ITEM}.priwig 
        mv ${IN}_${ITEM}.wig $TMPD/WIG/
    done
done

