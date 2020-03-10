#!/usr/bin/env bash
now=$(date +"%r")
echo "Current time : $now"

#############################
#sgRNA Designer txt to fasta#
#############################

#sgRNA_to_bed.sh Usage example:  sbatch --mem=32g -J sgBED --gres=lscratch:380  --time=30:00:00  --export=sgRNA_SET_NAME=SNAI2_MYOD_site4.2 sgRNA_to_bed.sh
# 
#options required:
# -sgRNA_SET_NAME 		Name of the file prefix used after downloading results from https://portals.broadinstitute.org/gpp/public/analysis-tools/sgrna-design

echo "> sgRNA set name: $sgRNA_SET_NAME"
sgRNAtxt=${sgRNA_SET_NAME}_sgrna-designs.txt
sgRNAstack=${sgRNA_SET_NAME}_sgrna.stack

echo "> forging $sgRNAstack"
cat $sgRNAtxt | cut -f 20 | sed -e '1d' > $sgRNAstack

echo "> looping BLAT for $sgRNAstack"
module load blat


for i in `cat $sgRNAstack`
do 
	if [ -f "$sgRNA_SET_NAME.$i.psl" ]
	then
		echo "> $sgRNA_SET_NAME.$i.psl blatted previously, skip to bed making."
	else
		echo "> making $sgRNA_SET_NAME.$i.fasta"
		##################
		#1-line fasta2bed# 
		##################
		printf ">sgRNA_$i_${sgRNA_SET_NAME}\n$i\n" >  $sgRNA_SET_NAME.$i.fasta
		
		echo "> converting to $sgRNA_SET_NAME.$i.nib"
		faToNib $sgRNA_SET_NAME.$i.fasta $sgRNA_SET_NAME.$i.nib

		echo "> blatting $sgRNA_SET_NAME.$i.nib into $sgRNA_SET_NAME.$i.psl"
		blat -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 /fdb/igenomes/Homo_sapiens/UCSC/hg19/hg19.fa $sgRNA_SET_NAME.$i.nib $sgRNA_SET_NAME.$i.psl
	fi
	
	echo "> cutting psl to make bed 4-column output!"
	cat $sgRNA_SET_NAME.$i.psl | sed '6q;d' | awk -v OFS='\t' '{print $14, $16, $17, $10, $9}' > $sgRNA_SET_NAME.$i.temp.bed
	strand=`cat $sgRNA_SET_NAME.$i.temp.bed | cut -f 5`
	echo "...sgRNA $i is on the $strand strand"
	
	if [ $strand = '+' ]
	then
		echo "...identifying cut site"
		cat $sgRNA_SET_NAME.$i.temp.bed |  awk -v OFS='\t' '{print $1, $2+16, $3-3, $4, $5}' | sed "s/.nib/.cut/" > $sgRNA_SET_NAME.$i.cutsite.bed
	fi
	if [ $strand = '-' ]
	then
		echo "...identifying cut site"
		cat $sgRNA_SET_NAME.$i.temp.bed |  awk -v OFS='\t' '{print $1, $2+3, $3-16, $4, $5}' | sed "s/.nib/.cut/" > $sgRNA_SET_NAME.$i.cutsite.bed
	fi	

done

echo "concatenating $sgRNA_SET_NAME.*.bed > $sgRNA_SET_NAME.bed"
cat $sgRNA_SET_NAME.*.temp.bed > $sgRNA_SET_NAME.sgRNA.bed
cat $sgRNA_SET_NAME.*.cutsite.bed > $sgRNA_SET_NAME.cutsites.bed

#sort -k1,1 -k2,2n $sgRNA_SET_NAME.bed > $sgRNA_SET_NAME.sort.bed
#cat $sgRNA_SET_NAME.sort.bed > $sgRNA_SET_NAME.bed

echo "> axing superfluous intermediate files"
rm -fv $sgRNA_SET_NAME.*.fasta $sgRNA_SET_NAME.*.nib $sgRNAstack $sgRNA_SET_NAME.*.temp.bed $sgRNA_SET_NAME.*.cutsite.bed $sgRNA_SET_NAME.*.psl

chgrp khanlab -R $sgRNA_SET_NAME*

later=$(date +"%r")
echo "Current time : $later"

echo "done, check results!"

#make a multi line fasta: cat $sgRNAtxt | cut -f 20 | sed 's/.*/&\n>guide/'| sed -e '1d' | sed '$d' > $sgRNAfasta
#faToTwoBit /fdb/fastadb/est_human.fas est_human.2bit
