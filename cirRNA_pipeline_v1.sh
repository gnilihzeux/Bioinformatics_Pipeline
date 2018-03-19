#!/bin/bash
# Copyright Statement
#
# Author
#   author: Jinwen
#   e-mail: zhangjinwen@gcbi.com.cn
# Description
#   aim:
#   inputs:
#   outputs:
# Software Dependencies
#   

# header ------------------------------------------------
d=`date +"%F-%T"`
exec >/data3/temp/$d.log 2>&1
workdir=$1
species=$2
group_file=$3

group_file=$workdir/$group_file

# softwares
PATH=/data2/linjie/rna_seq/soft/tophat-2.1.0/src/:$PATH
PATH=/opt/bowtie-1.0.0/:$PATH
PATH=/opt/bowtie2/:$PATH
PATH=/data2/linjie/circRNA/CIRCexplorer2-2.1.2/ucsc/:$PATH
PATH=/opt/software/:$PATH

# making dir 
if [[ ! ${workdir} =~ /$ ]]; then
  workdir=${workdir}/
fi

samples=`awk -F "\t" '{print $1}' $group_file`

for each in $samples;do
  mkdir -p -m 777 ${workdir}circRNA/${each}
done

for each in $samples;do
  mkdir -p -m 777 ${workdir}circRNA/${each}/circ_out/STAR
done

# main ---------------------------------------------------
for i in $samples;do
  
  STAR_dir=${workdir}circRNA/${i}/circ_out/STAR/
  FQ1=${workdir}Sample_${i}/${i}_combined_R1.fastq.gz
  FQ2=${workdir}Sample_${i}/${i}_combined_R2.fastq.gz
  DIR=${workdir}circRNA/${i}

  if [ "$species"=="human" ]; then
    # RUN STAR to align
    STAR --runThreadN 26 \
         --outFilterMultimapNmax 1 \
	     --genomeDir /data2/rna/star \
	     --outSAMtype BAM \
	     SortedByCoordinate \
	     --outFileNamePrefix $STAR_dir \
	     --chimSegmentMin 10 \
	     --outSAMstrandField intronMotif \
	     --readFilesCommand zcat \ 
	     --readFilesIn $FQ1 $FQ2

    cd $DIR
    # parse
    CIRCexplorer2M parse -t STAR ${workdir}circRNA/${i}/circ_out/STAR/Chimeric.out.junction > CIRCexplorer2_parse.log

    # assemble
    CIRCexplorer2M assemble -p 26 -r /data2/linjie/circRNA/CIRCexplorer2-2.1.2/test_4/hg38_ref.txt ${workdir}circRNA/${i}/circ_out > CIRCexplorer2_assemble.log

    # annotate
    CIRCexplorer2M annotate -r /data2/linjie/circRNA/CIRCexplorer2-2.1.2/test_4/hg38_ref.txt -g /data2/linjie/circRNA/CIRCexplorer2-2.1.2/test_4/hg38.fa ${workdir}circRNA/${i}/circ_out > CIRCexplor2_annotate.log

    # denovo
    CIRCexplorer2M denovo -r /data2/linjie/circRNA/CIRCexplorer2-2.1.2/test_4/hg38_ref.txt -g /data2/linjie/circRNA/CIRCexplorer2-2.1.2/test_4/hg38.fa ${workdir}circRNA/${i}/circ_out

    cd ..
  elif [ "$species"=="mmu" ]; then
    # RUN STAR to align
    STAR --runThreadN 15 \
	     --outFilterMultimapNmax 1 \
		 --genomeDir /data2/rna-mmu/star \
		 --outSAMtype BAM \
		 SortedByCoordinate 
		 --outFileNamePrefix $STAR_dir \
		 --chimSegmentMin 10 \
		 --outSAMstrandField intronMotif \
		 --readFilesCommand zcat \
		 --readFilesIn $FQ1 $FQ2

    cd $DIR
    # parse
    CIRCexplorer2M parse -t STAR circ_out/STAR/Chimeric.out.junction > CIRCexplorer2_parse.log

    # assemble
    CIRCexplorer2M assemble -p 15 -r /data3/project/circRNA_test/result/sample_21_2/mm10_ref.txt circ_out > CIRCexplorer2_assemble.log

    # annotate
    CIRCexplorer2M annotate -r /data3/project/circRNA_test/result/sample_21_2/mm10_ref.txt -g /data3/project/circRNA_test/result/sample_21_2/mm10.fa circ_out > CIRCexplor2_annotate.log

    # denovo
    CIRCexplorer2M denovo -r /data3/project/circRNA_test/result/sample_21_2/mm10_ref.txt -g /data3/project/circRNA_test/result/sample_21_2/mm10.fa circ_out

    cd ..
  fi

done

R --slave --vanilla --args $workdir $group_file hsa 1.2 1 0.5</data3/circRNA/bin/diff_circ_deseq_v1.R

