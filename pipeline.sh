#!/bin/bash
# NGS_PIPELINE
#This is a pipeline that analyzes chromosome 22 from a 3 year old boy with Myelodysplastic leukodystrophy (MLD).
#First we decompress 2 files containing paired and unpaired reads of chromosome 22 sequenced by HiSeq2500
#Assuming all conda tools are installed and directories are created:
#Go to the untrimmed_fastq directory to download files
echo the name of the script is $0
cd ~/ngs_course/dnaseq/data/untrimmed_fastq

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R1.fastq.gz

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R2.fastq.gz

#Download the bed file(contains genomic coordinates) in data
cd ..
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/chr22.genes.hg19.bed

#Then we need to download the Reference Human Genome file hg19.fa.gz
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

cd untrimmed_fastq
#fastq.gz is a compressed format. I need to uncompress it first
zcat WES01_chr22m_R1.fastq.gz > WES01_chr22m_R1.fastq
zcat WES01_chr22m_R2.fastq.gz > WES01_chr22m_R2.fastq

#In the searching_files lesson, some commands are covered like grep, less and sort to look for duplicate records and bad reads represented by 'NNNNNNNNNN' but we're going to use trimmomatic soon so I'm not including it here (Refer to workshop if needed)
#Now we're going to run fastqc on the raw untrimmed files
#Change directory to raw uncompressed fastq
#cd ~/ngs_course/dnaseq/data/untrimmed_fastq

fastqc -t 4 *.fastq.gz  #-t 4 is the multithreading functionality of fastqc to run faster

#move results to a new directory created within results
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

#Now take a closer look at the generated report from fastqc.
cd ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
ls -lh ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

#There are 2 outputs: one is a zip file and the other is an html report.
#Download Filezilla, setup using the username 'reem', ip address obtained via this command 'ip addr show eth0 | grep -oP '(?<=inet\s)\d+(\.\d+){3}' and the type of session as sftp (secure) and port 22. Click Quickconnect. The right hand side will be the wsl terminal while the left hand side will be the local PC. Choose a valid destination with easy permissions, eg. Downloads or Desktop. Transfer files by drag and drop or right-click and then download or view/edit. Both html reports are downloaded in Downloads. You can check the data.

#This is a for loop that iterates over the zip files from fastqc and unzips them
for zip in *.zip; do unzip $zip; done
ls -lh WES01_chr22m_R1.fastqc

#View the summary qc from both files (R1 and R2)
head WES01_chr22m_R1_fastqc/summary.txt
head WES01_chr22m_R2_fastqc/summary.txt

#concatenate both summary files into one file called fastqc_summaries and move it to logs
cat */summary.txt > ~/ngs_course/dnaseq/logs/fastqc_summaries.txt

#TRIMMOMATIC
#Run trimmomatic 0.39 to trim off bases that fall below a certain quality threshold, in this case 'phred33' using ILLUMINACLIP which cuts adapter and other illumina-specific sequences from the read.The TRAILING option cuts bases off the end of a read, if below a certain quality (25 in this case) and MINLEN drops an entire read if below a certain length (50bp). 
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
trimmomatic PE  \
 -threads 4 \
 -phred33 \
 /home/reem/ngs_course/dnaseq/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz /home/reem/ngs_course/dnaseq/data/untrimmed_fastq/WES01_> -baseout  /home/reem/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_trimmed_R \
 ILLUMINACLIP:/home/reem/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
 TRAILING:25 MINLEN:50

#After running trimmomatic we need to re-run fastqc on the trimmed files, so:
cd ~/ngs_course/dnaseq/data/trimmed_fastq
fastqc -t 4 *trimmed*

#Create a new diectory for the fastqc results of trimmed files and move them there
mkdir ~/ngs_course/dnaseq/results/fastqc_trimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_trimmed_reads/

#Then we move to the fastqc results folder we just created and unzip the compressed files
cd ~/ngs_course/dnaseq/results/fastqc_trimmed_reads
for zip in *.zip; do unzip $zip; done

#Concatenate summaries.txt in the same file as untrimmed files for review.
cat */summary.txt > ~/ngs_course/dnaseq/logs/fastqc_summaries.txt

#You can view the html reports on Filezilla as before and note the change in quality of reads after trimming.
#DONE!!!

#ALIGNMENT:
#This alignment is using bwa-mem.
#First, we build an index for the reference file to enable more efficient searching and save space.

cd ~/ngs_course/dnaseq/data/reference
#Because I only have 3.8 Gb of space on WSL I could not index the whole reference genome hg19, so I downloaded the reference chromosome 22 only from the hg19 build.
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz

#Now it's time to index the ref. chr22
#bwa index ~/ngs_course/dnaseq/data/reference/chr22.fa.gz

#Create a directory for all the aligned data to be processed as follows:
mkdir ~/ngs_course/dnaseq/data/aligned_data

#Run bwa-mem using RG Info. Read Group information provides essential metadata that guides the alignment process
#and ensures that sequencing data is correctly attributed to its source.
#This metadata is essential for downstream analyses and maintaining data quality in genomics research.
#The output is a 'sam' file.
#file to save space and be in a more reader-friendly format, the bam file is then sorted by reference sequence and position of aligned reads.
#This allows for efficient retrieval of reads corresponding to specific genomic regions.
#An index file (.bam.bai) allows for efficient retrieval of reads corresponding to specific genomic regions.
cd ~/ngs_course/dnaseq/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/chr22.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m.sam
samtools view -h -b WES01_chr22m.sam > WES01_chr22m.bam

samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam

samtools index WES01_chr22m_sorted.bam #This will generate a .bai index file

#Check the outputs. You should have a sam, bam, sorted bam as well as indexed bam.
ls

#POST-ALIGNMENT QC AND FILTERING:
#Mark duplicate reads using picard, you can view the txt file using less or download it through Filezilla
picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt

#The marked file is then indexed using samtools as before
samtools index WES01_chr22m_sorted_marked.bam

#Next:
#samtools view: This is a command from the samtools toolkit for processing SAM/BAM files.
#-F 1796: a filter that excludes reads with specific flag bits set. In this case, 1796 is a decimal number, but it's often specified in binary. In binary, 1796 is 11100000100. This indicates that you want to filter out reads with the following flag bits set:
#In short, you are excluding reads that are unmapped, unmapped as part of a pair, fail quality checks, or are supplementary alignments.
samtools view -F 1796  -q 20 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam

#The sorted, filtered, marked bam is then indexed as before.
samtools index WES01_chr22m_sorted_filtered.bam


#Samtools flagstat provides statistics about the alignment flags and mapping of reads in a BAM file (whether it's paired, properly mapped or secondary alignment) to assess quality and characteristics before downstream analysis.
samtools flagstat WES01_chr22m_sorted_filtered.bam

#Samtools idxstat provides statistics about the alignment of reads mapped to different reference sequences. It is useful to know how many reads are mapped to specific regions of the reference genome.
samtools idxstats WES01_chr22m_sorted_filtered.bam

#Picard calculates the length of the DNA fragment being sequenced and computes various metrics (eg. mean, median, standard deviation) as well as generates a histogram with the most common insert sizes as a pdf file.
picard CollectInsertSizeMetrics \
I=WES01_chr22m_sorted_filtered.bam \
O=insert_size_metrics.txt \
H=insert_size_histogram.pdf \
M=0.5

#Bedtools calculates the depth of coverage which is the number of bases covered by reads in a specific region. This helps identify regions with low coverage.
bedtools coverage -a /home/reem/ngs_course/dnaseq/data/chr22.genes.hg19.bed -b WES01_chr22m_sorted_filtered.bam > coverage_output.txt

#You can view your sorted and filtered BAM file using IGV by uploading both the '.bam' and '.bam.bai' file. Zoom in to see reads with more coverage vs low coverage as well as identify regions where there are indels or substitutions. 
#DONE

#VARIANT CALLING WITH FREEBAYES:
#First we unzip the reference genome
#cd ~/ngs_course/dnaseq/data/reference
#zcat chr22.fa.gz > chr22.fa

#Then, samtools faidx is used to index the referece genome in FASTA format (for faster retrieval of sequences as mentioned before)
samtools faidx chr22.fa

#Run 'Freebayes' variant caller to call variants from the aligned bam file using the reference genome chr22.fa to which the reads will be compared to call variants. The output is a vcf file.
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/chr22.fa --vcf ~/ngs_course/dnaseq/results/WES01_chr22m.vcf


#The vcf file is then compressed using bgzip for more efficient storage and transfer. (Tabix has to be installed before this step)
cd ../..
cd results
bgzip ~/ngs_course/dnaseq/results/WES01_chr22m.vcf

#Finally, tabix indexes the compressed VCF file for subsequent analysis and querying
tabix -p vcf ~/ngs_course/dnaseq/results/WES01_chr22m.vcf.gz

#FILTERING:
#vcflib is a tool that filters called variants according to parameters specified in the command:
#exclude low quality calls (allow QUAL>1), quality to alternate allele ratio < 10, low read depth, calls where observations are observed on only one strand (forward (SAF) or reverse(SAR)), unbalanced reads where observations are made only when reads are placed either to the left(5') or right(3') of the alternate allele (RPL>1 and RPR>1)
#conda install vcflib #if you did not install vcflib before. vcffilter is part of the vcflib suite

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq/results/WES01_chr22m_filtered.vcf

#Bedtools intersect performs an intersection operation between 2 files: the filtered vcf file and the bed file which contains the genes targeted in this trial data. -header: includes the header from the input vcf file in the output. -wa: writes the original entries in file A( the vcf file) that overlap with file B (the bed file).
#This results in an output vcf file (WES01_chr22m_filtered_chr22.vcf) which is redirected to the results directory
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/WES01_chr22m_filtered.vcf -b /home/reem/ngs_course/dnaseq/data/chr22.genes.hg19.bed   > ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf

#The following commands just compress the resulted output vcf file (.gz) and indexes it (.tbi) using tabix as before.
bgzip ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf

tabix -p vcf ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf.gz

#ANNOTATION:
#This step involves using annovar, a software for annotating called variants and generating a csv file.
#Download annovar by filling the registration form to get a download link
#Extract the tar.gz file downloaded
#tar -zxvf annovar.latest.tar.gz    #If not already done
#Download annovar databases
#cd annovar
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/


#Assuming the previous steps have already been done, we will proceed to converting the vcf file to an annovar compatible file:
#Now, move to the directory where annovar was installed
cd /mnt/c/Users/MOSTAFA/Downloads/annovar

#This command takes a vcf file version 4 as input and converts it to an avinput file using convert2annovar.pl script and saves it in the specified directory
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf.gz > ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.avinput


#./table_annovar is the script we're running on the WES01_chr22m_filtered_chr22.avinput query file using the database humandb generating an output 'WES01_chr22m_filtered_chr22' and -remove just tells annovar to remove temporary files after annotation. -protocol specifies the annotation protocols used, i.e. refGene, ensGene, clinvar....etc.. Operation g stands for gene-based annotations while f stands for filter-based annotations. -nastring specifies the string to use when a score is not available (here it's empty so don't put any value) and finally -csvout means generate a comma_delimited csv output file.
./table_annovar.pl ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19 \
   -out ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22 -remove \
   -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

#Install snpeff: I'm going to assume this is already done!
#cd
#conda install snpeff

#To search for snpeff after installing to locate the jar file:
#find / -name snpEff
#My snpeff location:$ cd /home/reem/anaconda3/bin/

#Install the hg19 database (GRCh37.75)
#java -jar /home/reem/anaconda3/pkgs/snpeff-5.1-hdfd78af_0/share/snpeff-5.1-0/snpEff.jar download -v GRCh37.75
#ls
#mv ~/anaconda3/pkgs/snpeff-5.1-hdfd78af_0/share/snpeff-5.1-0/data/GRCh37.75 ~/ngs_course/dnaseq/data/reference/
#I realized you cannot conduct snpeff using a vcf file that was aligned against chr22.fa while snpeff will use hg19..........Different reference genome databases.Therefore, I repeated the alignment with hg19 again, since I managed to allocate more memory to wsl. (8.4Gb instead of 3.7Gb)
#The following commands repeat the exact same steps as with chr22.fa but with hg19.fa.gz
#cd ngs_course/dnaseq/data/
#cd reference/
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa
cd ~/ngs_course/dnaseq/data/aligned_data
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/WES01_hg19.vcf
tabix -p vcf ~/ngs_course/dnaseq/results/WES01_hg19.vcf.gz
ls
cd ../..
cd results/
bgzip ~/ngs_course/dnaseq/results/WES01_hg19.vcf
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq/results/WES01_hg19.vcf.gz > ~/ngs_course/dnaseq/results/WES01_hg19_filtered.vcf

#Now we can conduct snpeff analysis using the new vcf 'WES01_hg19.vcf'. At this point, I ran out of memory so I deleted the raw fastq.gz files.
cd ..
cd data/untrimmed_fastq
rm -r *.gz
cd ..

Create a new directory to conduct the snpeff analysis. -Xmx4G in the sneff command helps allocate more memory to the JVM (Java Virtual Machine), 4G in this case to help run snpeff. 
mkdir snpeff
cd snpeff/
java -Xmx4G -jar /home/reem/anaconda3/pkgs/snpeff-5.1-hdfd78af_0/share/snpeff-5.1-0/snpEff.jar eff -dataDir /home/reem/ngs_course/dnaseq/data/reference/ -v GRCh37.75 /home/reem/ngs_course/dnaseq/results/WES01_hg19_filtered.vcf > snpeff_annotations.vcf

#You can view the snpeff_annotations.vcf file on IGV. Now we're going to filter to exonic variants only.
grep 'ANN=.*protein_coding.*' snpeff_annotations.vcf > exonic_variants.vcf

#If you want to check the number of variants in the final output, you can use this command
grep -E 'ANN=.*protein_coding.*' exonic_variants.vcf | wc -l

#In this step, I filter the annovar output to produce only exonic variants not seen in ExAC (Exome Aggregation Consortium). (It was supposed to be dbSNP but there's no dbSNP in this dataset).
cd ../..
cd results/
#cat WES01_hg19_filtered.vcf
awk -F, 'NR == 1 || ($21 == "." && $22 == "." && $23 == "." && $24 == "." && $25 == "." && $26 == "." && $27 == "." && $28 == "." && $6 ~ /exonic/) { print }' WES01_chr22m_filtered_chr22.hg19_multianno.csv > exonic_variants_not_in_Exac.csv

#The previous command looks complicated but it actually isn't. 'awk' is a command used to filter csv files. '-F,' specifies that the field separator as comma ','. 'NR == 1' is the Number of Row == 1, i.e. if it's the first (containing the headers), then yes it's included in the filtering criteria. '||' is the OR Operator specifying that if both conditions on either side of the '||' are true, then proceed to the next part of the command. ($21 == "." && $22 == "." && $23 == "." && $24 == "." && $25 == "." && $26 == "." && $27 == "." && $28 == "." This part corresponds to the columns that match to ExAC in the csv file. If these columns are empty, indicated by a period ('.'), meaning that the variants have not been seen in ExAC plus $6, which corresponds to the functional class of the gene, is an exon, then {print} and redirect (>) the output to a new file 'exonic_variants_not_in_Exac.csv'.

#THE END
#REEM HASSAN






















