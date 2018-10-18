## declare an array variable
declare -a arr=("IIIweekCtrl_S4_R1_001-TR.fastq.gz" "IIIweekEV_S5_R1_001-TR.fastq.gz" "IIIweekRDR_S6_R1_001-TR.fastq.gz" "IIweekCtrl_S1_R1_001-TR.fastq.gz" "IIweekEV_S2_R1_001-TR.fastq.gz" "IIweekRDR_S3_R1_001-TR.fastq.gz")
# now loop through the above array
for file in "${arr[@]}"
do
	echo "$file"
	gunzip $file
done
