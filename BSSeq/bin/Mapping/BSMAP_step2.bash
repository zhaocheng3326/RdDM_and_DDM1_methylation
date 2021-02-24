#!/bin/bash
DIR=/cluster/home/lihe/My_project/BSSeq
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin


Rscript $SRC/methy_barplot_report.R $TMPD $RE/meth_summary.temp Total_region
Rscript $SRC/Gene_methylation_PCA.R

/psc/program/install/R-3.1.1/bin/Rscript $SRC/methy_barplot_report.R $TMPD/  $RE/meth_summary.temp Total_region
/psc/program/install/R-3.1.1/bin/Rscript $SRC/methy_barplot_report.R $TMPD/  $TMPD/IG.out.summary.temp  Intergenic
/psc/program/install/R-3.1.1/bin/Rscript $SRC/methy_barplot_report.R $TMPD/  $TMPD/protein_coding_gene_ML.out.summary.temp  Protein_coding
/psc/program/install/R-3.1.1/bin/Rscript $SRC/methy_barplot_report.R $TMPD/  $TMPD/TE.out.summary.temp TE


