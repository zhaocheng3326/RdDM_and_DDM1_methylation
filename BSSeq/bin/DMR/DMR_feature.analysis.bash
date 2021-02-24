#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq/
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin





cd $TMPD
echo -e "ID\tLocation\tType\tCol\tpol4\tpol5\tddm1\td4\td5\tcmt2\tddm1cmt2" > $TMPD/TE.nearBy.60bin.out
paste <(paste <(cat Col_1/Col_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat Col_2/Col_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat pol4_1/pol4_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat pol4_2/pol4_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat pol5_1/pol5_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat pol5_2/pol5_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat ddm1_1/ddm1_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat ddm1_2/ddm1_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat d4_1/d4_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat d4_2/d4_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat d5_1/d5_1_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat d5_2/d5_2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}')  <(cat cmt2/cmt2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat ddm1cmt2/ddm1cmt2_TE.nearBy.60bin.out |perl $SRC/cal_methy.pl) |cut -f 1,2,3,6,9,12,15,18,21,24 |sed -e 's/_/\t/g'|cut -f 1,3,4,5,6,7,8,9,10,11,12 >> $TMPD/TE.nearBy.60bin.out


echo -e "ID\tLocation\tType\tCol\tpol4\tpol5\tddm1\td4\td5\tcmt2\tddm1cmt2" > $TMPD/Gene.nearBy.60bin.out
paste <(paste <(cat Col_1/Col_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat Col_2/Col_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat pol4_1/pol4_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat pol4_2/pol4_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat pol5_1/pol5_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat pol5_2/pol5_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat ddm1_1/ddm1_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat ddm1_2/ddm1_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat d4_1/d4_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat d4_2/d4_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(paste <(cat d5_1/d5_1_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat d5_2/d5_2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl)|awk '{OFS="\t";print $1,$2,($3+$6)/2}') <(cat cmt2/cmt2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) <(cat ddm1cmt2/ddm1cmt2_gene.nearBy.60bin.out |perl $SRC/cal_methy.pl) |cut -f 1,2,3,6,9,12,15,18,21,24 |sed -e 's/_/\t/g'|cut -f 1,3,4,5,6,7,8,9,10,11,12 >> $TMPD/Gene.nearBy.60bin.out



DMR=$TMPD/DMR_CHH_logFC_4_1.58


cd $DMR
#### building exact region


#### make control region
shuffleBed -i total.analysis.QS.bed -chrom -seed 1  -g ~/Genome/ATH/annotation/TAIR10.bg -excl total.analysis.QS.bed|awk -F "\t" '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3"_control"}' > control.region.bed


#### The size of overlaped TE
for ITEM in C1 C2 C3  control.region ; do  intersectBed -a ${ITEM}.bed  -b  <(cat $DOC/TE_bedlike.bed|awk '{if (($3-$2) >=100) print $0}') -wa -wb |cut -f 4,8,9|sort|uniq  > ${ITEM}.OV.TE.size ; done


#### miRNA target
for ITEM in C1 C2 C3  control.region ; do intersectBed -a ${ITEM}.bed  -b <(cat $DOC/table.S3.TE_miRNA.out |cut -f 1-4|cut -f 1 -d "_"|sed -e 's/ath-//') -wa -wb |cut -f 1-8 |sort |uniq > ${ITEM}.OV.miRNA;done

#### DMR coverage during the 60 bin of gene and te
mkdir $DMR/Asso_Gene/
mkdir $DMR/Asso_TE/
for ITEM in C1 C2 C3 control.region
do
    cd $DMR

	intersectBed -a <(awk -F "\t" '{OFS="\t";print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TAIR10.C.gene.bed|perl $SRC/UD_modify.pl) -b  <(awk -F "\t" '{OFS="\t";print $1,$2,$3,$4,".\t+"}' ${ITEM}.bed) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g'|awk -F "\t" '{OFS="\t";print $2"|"$3"|"$4,$1,$6}' |sort|uniq  > Asso_Gene/${ITEM}.Inter.Gene.out

	intersectBed -a <(awk -F "\t" '{OFS="\t";if  (($3-$2) >99) print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TE_bedlike.bed|perl $SRC/UD_modify.pl) -b  <(awk -F "\t" '{OFS="\t";print $1,$2,$3,$4,".\t+"}' ${ITEM}.bed) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g'|awk -F "\t" '{OFS="\t";print $2"|"$3"|"$4,$1,$6}' |sort|uniq  > Asso_TE/${ITEM}.Inter.TE.out
done


#### Annotation of DMR
for ITEM in C1 C2 C3 control.region
do
    cd $DMR
    cat <(intersectBed  -a ${ITEM}.bed -b $DOC/Protein_coding_gene.bed -wa -wb |cut -f 4 |sort|uniq| awk '{print $0"\tPC"}') <(intersectBed  -a ${ITEM}.bed -b <(cat $DOC/TE_bedlike.bed|awk '{if (($3-$2) >99) print $0}') -wa -wb |cut -f 4 |sort|uniq| awk '{print $0"\tTE"}') |sort -k1,1 |groupBy -g 1 -c 2 -o collapse |sed -e 's/,/_/' |perl $SRC/buqi_region.pl ${ITEM}.bed|awk -v a=${ITEM} '{print $0"\t"a}' > ${ITEM}.type
done

cat C1.type C2.type C3.type control.region.type  >  class_total.type


#### Annotation of DMR (bp based)
for ITEM in C1 C2 C3 control.region
do
    cd $DMR
    mkdir -p BP_temp
    intersectBed  -a  <(bedtools makewindows -b ${ITEM}.bed -w 1|sort |uniq) -b <(cat <(cat $DOC/TE_bedlike.bed|awk '{OFS="\t";print $1,$2,$3,"4"}') <(cat $DOC/Protein_coding_gene.bed|awk '{OFS="\t";print $1,$2,$3,"3"}')  <(cat $DOC/Protein_coding_gene.upstream.bed|awk '{OFS="\t";print $1,$2,$3,"2"}') ) -wa -wb |cut -f1-3,7|sort -k1,1 -k2,2 -k3,3 |groupBy -g 1,2,3, -c 4  -o max > BP_temp/${ITEM}.split.bp.bed
    echo -e "Total_bp\t`bedtools makewindows -b ${ITEM}.bed -w 1|awk '{print $3-$2}'|SM`" > ${ITEM}.bp.anno
    cut -f 4 BP_temp/${ITEM}.split.bp.bed|perl ~/PC/code/count_repeat.pl |sed -e 's/^2/Promoter/'|sed -e 's/^3/PC_gene/'|sed -e 's/^4/TE/' >> ${ITEM}.bp.anno
done
\rm -r $DMR/BP_temp




#### The distribution of DMR on whole chromosome
for ITEM in C1 C2 C3 control.region
do
	intersectBed -a $DOC/TAIR10.bg.100kp.windows.bed -b ${ITEM}.bed -wa -c |grep -v "ChrC"|grep -v "ChrM" |cut -f 4,5|sed -e 's/_/\t/' > ${ITEM}.DMR.100kb.cov
done



#### GC information of DMR
# GC content CpG content
# SF feature GC content=(N_C+N_G)/L ,  CpG 0/E  = N_CpG*L/(Nc*NG) GC and CpG

for ITEM in C1 C2 C3 control.region; do  join -1 4 -2 1 -t$'\t' <(bedtools nuc -fi $DOC/arab10.C.fa -bed ${ITEM}.bed  |grep -v "^#"|cut -f 1-4,6,8,9 |sort -k4,4)  <(intersectBed -a ${ITEM}.bed   -b <(grep "CG" ~/Genome/ATH/DNA/BG_forw.txt|awk '{OFS="\t";print $1,$2,$2+1,$4}') -c |cut -f4,5|sort -k1,1)|awk '{OFS="\t";if ($6*$7 !=0) print $1,$5,$8*($4-$3)/($6*$7)}' > ${ITEM}.GC.contend ; done


# Distance between DMR and nearest gene TSS
for ITEM in C1 C2 C3  control.region; do closestBed -a <(cat ${ITEM}.bed |sortBed) -b <(cat $DOC/Protein_coding_gene.bed|awk '{OFS="\t";if ($6=="+") print $1,$2,$2+1,$4,"TSS",$6;else print $1,$3-1,$3,$4,"TSS",$6}'|sortBed) -t first -d |cut -f 4,8,11   > ${ITEM}.closed_TSS.dis ; done

# Fraction of DMR overlapping with gene promoters ##ID 1K 2K
for ITEM in C1 C2 C3 control.region; do join -1 1 -2 1 -t$'\t' <(intersectBed -a ${ITEM}.bed  -b <(cut -f1-6 $DOC/Protein_coding_gene.bed |awk '{OFS="\t";if ($6=="+") print $1,$2-1000,$2,$4,"Upstream",$6;else print $1,$3,$3+1000,$4,"Upstream",$6}') -wa -c |awk '{if ($5 =="0") print $4"\tNon_overlap";else print $4"\tOverlap"}' |sort -k1,1)  <(intersectBed -a ${ITEM}.bed  -b $DOC/Protein_coding_gene.upstream.bed -wa -c |awk '{if ($5 =="0") print $4"\tNon_overlap";else print $4"\tOverlap"}'|sort -k1,1)  > ${ITEM}.promoter.1k_2k.ov; done

# DMR over TEs of different sizes #region_code(TE 1 means -4000+100*1)
for ITEM in C1 C2 C3 control.region;  do intersectBed -a <(join -1 1 -2 4 -t$'\t' <(intersectBed -a <(cat $DOC/TE_bedlike.TSS_nearBy4k.100inter.bed|sed -e 's/_/\t/g') -b ${ITEM}.bed  -wa -wb |cut -f 4 |sort |uniq |sort -k1,1 ) <(cat $DOC/TE_bedlike.TSS_nearBy4k.100inter.bed|sed -e 's/_/\t/g'|sort -k4,4 ) |awk '{OFS="\t";print $2,$3,$4,$1"_"$6}') -b ${ITEM}.bed  -wa -wb|cut -f 4,8|sed -e 's/_/\t/' > ${ITEM}.TSS_nearBy.cov; done
for ITEM in C1 C2 C3 control.region;  do intersectBed -a <(join -1 1 -2 4 -t$'\t' <(intersectBed -a <(cat $DOC/TE_bedlike.TTS_nearBy4k.100inter.bed|sed -e 's/_/\t/g') -b ${ITEM}.bed  -wa -wb |cut -f 4 |sort |uniq |sort -k1,1 ) <(cat $DOC/TE_bedlike.TTS_nearBy4k.100inter.bed|sed -e 's/_/\t/g'|sort -k4,4 ) |awk '{OFS="\t";print $2,$3,$4,$1"_"$6}') -b ${ITEM}.bed  -wa -wb |cut -f 4,8|sed -e 's/_/\t/' > ${ITEM}.TTS_nearBy.cov; done

# 20+20+20 (2kb+TE_body+2kb) DMR density
for ITEM in C1 C2 C3 control.region;  do intersectBed -a <(awk -F "\t" '{OFS="\t";if  (($3-$2) >99) print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TE_bedlike.bed|perl $SRC/UD_modify.pl) -b  <(awk -F "\t" '{OFS="\t";print $1,$2,$3,$4,".\t+"}' ${ITEM}.bed) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g'|awk -F "\t" '{OFS="\t";print $1,$6,$2"|"$3"|"$4}' |sort|uniq > ${ITEM}.60bin.cov; done

# pol4 coverage
#for ITEM in C1 C2 C3 control.region;  do intersectBed -a ${ITEM}.bed -b $DATA/pol4_chip/pol4mt.log2output.bed -wa -wb |cut -f 1-4,8-10 > ${ITEM}.pol4.50bp.out;done
#for ITEM in C1 C2 C3 control.region;  do intersectBed -a ${ITEM}.bed -b $DATA/pol5_chip/pol5mt.log2output.bed -wa -wb |cut -f 1-4,8-10 > ${ITEM}.pol5.50bp.out;done

# 20+20+20 (2kb+TE_body+2kb) Methylation_level density
for ITEM in C1 C2 C3 control.region;do head -n1  $TMPD/TE.nearBy.60bin.out |cut -f 2-11 > ${ITEM}.60bin.ML.out; join -1 1 -2 1 -t$'\t' <(cut -f 2 Asso_TE/${ITEM}.Inter.TE.out|sort -k1,1) <(cat $TMPD/TE.nearBy.60bin.out|sort -k1,1 )|sort -k2,2 -k3,3|groupBy  -g 2,3 -c 4,5,6,7,8,9,10,11 -o mean,mean,mean,mean,mean,mean,mean,mean >> ${ITEM}.60bin.ML.out;done

for ITEM in C1 C2 C3 control.region;do head -n1  $TMPD/TE.nearBy.60bin.out |cut -f 2-11|awk '{print "LenC\t"$0}' > ${ITEM}.60bin.sep.ML.out; join -1 1 -2 1 -t$'\t' <(join -1 1 -2 1 -t$'\t' <(cut -f 2 Asso_TE/${ITEM}.Inter.TE.out|sort -k1,1 ) <(cat $DOC/TE_length.cut.out|awk '{if ($2 >=100) print $0}'|cut -f 1,3|sort -k1,1)|sort -k1,1) <(cat $TMPD/TE.nearBy.60bin.out|sort -k1,1 )|sort -k2,2 -k3,3 -k4,4|groupBy  -g 2,3,4 -c 5,6,7,8,9,10,11,12 -o mean,mean,mean,mean,mean,mean,mean,mean >> ${ITEM}.60bin.sep.ML.out;done

join -1 1 -2 1  -t $'\t' <(intersectBed -a <(awk -F "\t" '{OFS="\t";print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TAIR10.C.gene.bed|perl $SRC/UD_modify.pl) -b  <(cat C1.bed C2.bed  C3.bed control.region.bed |awk -F "\t" '{OFS="\t";print $1,int(($2+$3)/2),int(($2+$3)/2)+1,$4,".\t+"}' ) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g' |awk '{OFS="\t";print $1,$6,$2"|"$3"|"$4}'|sort -k1,1)  <(sort -k1,1 $DOC/Protein_coding_gene.ID) > DMR.summit.ov.PC.gene.out

intersectBed -a <(awk -F "\t" '{OFS="\t";if  (($3-$2) >99) print $1,$2-2000,$3+2000,$4,$5,$6}' $DOC/TE_bedlike.bed|perl $SRC/UD_modify.pl) -b  <(cat C1.bed C2.bed C3.bed control.region.bed |awk -F "\t" '{OFS="\t";print $1,int(($2+$3)/2),int(($2+$3)/2)+1,$4,".\t+"}' ) -wa -wb |perl $SRC/cov_cal_UD.pl 2000 60 |sed -e 's/|/\t/g' |awk '{OFS="\t";print $1,$6,$2"|"$3"|"$4}'  > DMR.summit.ov.TE.out

### check the DMR distribution on the edge of TE
#cd $DOC
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-400,$2+400,$4,"leftE",$6,$7,$8"\n"$1,$3-400,$3+400,$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+400,$3-400,$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed

#cd $DOC
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-int($5*0.1),$2+int($5*0.1),$4,"leftE",$6,$7,$8"\n"$1,$3-int($5*0.1),$3+int($5*0.1),$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+int($5*0.1),$3-int($5*0.1),$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed
#cd $DOC
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-500,$2+500,$4,"leftE",$6,$7,$8"\n"$1,$3-500,$3+500,$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+int($5*0.25),$3-int($5*0.25),$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed

#cd $DOC
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-int($5*0.25),$2+int($5*0.25),$4,"leftE",$6,$7,$8"\n"$1,$3-int($5*0.25),$3+int($5*0.25),$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+int($5*0.25),$3-int($5*0.25),$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed

#cd $DOC
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-int($5*0.2),$2+int($5*0.2),$4,"leftE",$6,$7,$8"\n"$1,$3-int($5*0.2),$3+int($5*0.2),$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
#cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+int($5*0.2),$3-int($5*0.2),$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed

# not sure
cd $DOC
cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2-int($5*0.15),$2+int($5*0.15),$4,"leftE",$6,$7,$8"\n"$1,$3-int($5*0.15),$3+int($5*0.15),$4,"rightE",$6,$7,$8}' > m4k.TE_bedlike.edgeBd.bed
cat TE_bedlike.bed |awk '{if ($5>4000) print $0}'|awk -F "\t" '{OFS="\t";print $1,$2+int($5*0.15),$3-int($5*0.15),$4,"Body",$6,$7,$8}' >> m4k.TE_bedlike.edgeBd.bed

cd ~/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58
intersectBed -a  total.analysis.QS.bed -b $DOC/m4k.TE_bedlike.edgeBd.bed -wa -wb |cut -f 4,8,9|sort|uniq |sed -e 's/|/\t/g'|cut -f 4,6 |sort |uniq -c

# DMR distribution on TE (more than 4k)
for ITEM in C1 C2 C3 control.region;do  intersectBed -a  ${ITEM}.bed -b $DOC/m4k.TE_bedlike.edgeBd.bed -wa -wb |cut -f 4,8,9|sort|uniq > ${ITEM}.4kTE.bodyedge.out;done
#
#Rscript $SRC/DMR_feature_analysis.new.R









### calculate some extra region from the overlap.tiff
cd /cluster/home/lihe/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58/tempB

awk -F "\t"  -v a="(cmt2∩d45)∖(ddm1∪p45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  > temp.cmt2.merge.bed
awk -F "\t"  -v a="(cmt2∩ddm1∩d45)∖(p45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  >> temp.cmt2.merge.bed
awk -F "\t"  -v a="(cmt2)∖(ddm1∪p45∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  >> temp.cmt2.merge.bed
awk -F "\t"  -v a="(cmt2∩ddm1)∖(p45∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  >> temp.cmt2.merge.bed
awk -F "\t"  -v a="(cmt2∩p45∩d45)∖(ddm1)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  >> temp.cmt2.merge.bed
awk -F "\t"  -v a="(cmt2∩p45)∖(ddm1∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2"}' VD.output.merge.bed  >> temp.cmt2.merge.bed

awk -F "\t"  -v a="(cmt2∩p45∩d45)∖(ddm1)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  > temp.p45.merge.bed
awk -F "\t"  -v a="(cmt2∩p45)∖(ddm1∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  >> temp.p45.merge.bed
awk -F "\t"  -v a="(ddm1∩p45)∖(cmt2∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  >> temp.p45.merge.bed
awk -F "\t"  -v a="(ddm1∩p45∩d45)∖(cmt2)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  >> temp.p45.merge.bed
awk -F "\t"  -v a="(p45)∖(cmt2∪ddm1∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  >> temp.p45.merge.bed
awk -F "\t"  -v a="(p45∩d45)∖(cmt2∪ddm1)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"p45"}' VD.output.merge.bed  >> temp.p45.merge.bed

awk -F "\t"  -v a="(cmt2∩p45∩d45)∖(ddm1)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2p45"}' VD.output.merge.bed > temp.cmt2p45.merge.bed
awk -F "\t"  -v a="(cmt2∩p45)∖(ddm1∪d45)" '{OFS="\t"; if ($5 ==a) print $1,$2,$3,"cmt2p45"}' VD.output.merge.bed >> temp.cmt2p45.merge.bed


echo -n "" > ../extra.total.analysis.200version.bed
for ITEM in cmt2 p45 cmt2p45
do
        cd /cluster/home/lihe/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58/tempB
        sortBed -i temp.${ITEM}.merge.bed |mergeBed|awk -F "\t" -v a=${ITEM} '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3"\t"a}' >> ../extra.total.analysis.200version.bed
done

cd /cluster/home/lihe/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58
intersectBed -a <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 $TMPD/Col_meCHH.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk '{OFS="\t";print $4,$2,$3}') -b extra.total.analysis.200version.bed -f 0.99 -wa -wb |sort -k7,7 -k8,8 -k1,1 |groupBy -g 1,7,8 -c 2,3 -o min,max |awk '{OFS="\t";print $1,$4,$5,$1"|"$4"|"$5"|"$3,$3}'  > extra.total.analysis.QS.temp

echo -n "" > extra.total.analysis.QS.bed
for CLA in cmt2 p45 cmt2p45
do
        awk -F "\t" -v a=${CLA} '{OFS="\t";if ($5==a) print $1,$2,$3,$4}' extra.total.analysis.QS.temp |sortBed |mergeBed -d 50 |awk -v a=${CLA}  '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3"|"a}' |awk '{if (($3-$2)>= 50) print $0}' >> extra.total.analysis.QS.bed
done

cut -f 4 extra.total.analysis.QS.bed|cut -f 4 -d "|"|sort |uniq -c



#### not sure whether used or not
## Don't consider any overlap from the tiff file, use the hypo DMR directly
cd ~/My_project/BSSeq/tmp_data/DMR_CHH_logFC_4_1.58/tempB
intersectBed -a <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 $TMPD/Col_meCHH.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk '{OFS="\t";print $4,$2,$3}') -b <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 p45.hypo.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk -F "\t" '{OFS="\t";print $3,$4,$5,$1}'|sortBed |mergeBed |awk -F "\t" '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3}' ) -f 0.99 -wa -wb|sort -k7,7 -k1,1 |groupBy -g 1,7 -c 2,3 -o min,max |awk '{OFS="\t";print $1,$3,$4,$1"|"$3"|"$4}'|awk '{if (($3-$2)>= 50) print $0}'  > p45.hypo.out.QS.bed

intersectBed -a <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 $TMPD/Col_meCHH.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk '{OFS="\t";print $4,$2,$3}') -b <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 cmt2.hypo.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk -F "\t" '{OFS="\t";print $3,$4,$5,$1}'|sortBed |mergeBed |awk -F "\t" '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3}' ) -f 0.99 -wa -wb|sort -k7,7 -k1,1 |groupBy -g 1,7 -c 2,3 -o min,max |awk '{OFS="\t";print $1,$3,$4,$1"|"$3"|"$4}' |awk '{if (($3-$2)>= 50) print $0}' > cmt2.hypo.out.QS.bed

intersectBed -a <(join -1 1 -2 4 -t$'\t' <(sort -k1,1 $TMPD/Col_meCHH.out) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk '{OFS="\t";print $4,$2,$3}') -b <(join -1 1 -2 4 -t$'\t' <(cat cmt2.hypo.out p45.hypo.out |sort|uniq -d |sort -k1,1) <(sort -k4,4 $DOC/TAIR10.windows.bed )|awk -F "\t" '{OFS="\t";print $3,$4,$5,$1}'|sortBed |mergeBed |awk -F "\t" '{OFS="\t";print $1,$2,$3,$1"|"$2"|"$3}' ) -f 0.99 -wa -wb|sort -k7,7 -k1,1 |groupBy -g 1,7 -c 2,3 -o min,max |awk '{OFS="\t";print $1,$3,$4,$1"|"$3"|"$4}' |awk '{if (($3-$2)>= 50) print $0}' > cmt2p45.hypo.out.QS.bed
