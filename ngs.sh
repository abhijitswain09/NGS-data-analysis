#Quality analysis of Raw Data
fastqc */*.fastq.gz
#Merge all the qc reports
multiqc .
#Trim the Raw Data
#trim.sh = java -jar /var/abhi/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 6 /var/abhi/$1/$1\_1.fastq.gz /var/abhi/$1/$1\_2.fastq.gz /var/abhi/$1/$1\_output_forward_paired.fq.gz /var/abhi/$1/$1\_output_forward_unpaired.fq.gz /var/abhi/$1/$1\_output_reverse_paired.fq.gz /var/abhi/$1/$1\_output_reverse_unpaired.fq.gz ILLUMINACLIP:/var/abhi/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
ls -d SRR107526*|sort |while read i;do sh trim.sh $i ;done
#QC of the output files
fastqc */*.fq.gz
#Merge the QC report
multiqc .
#Build the Index File of the Genome
hisat2-build -f GCF_004153795.1_AHAU_CSS_1_genomic.fna Tea.index
gunzip */*_*_forward_paired.fq.gz */*_*_reverse_paired.fq.gz
ls -d SRR107526*|sort|while read i;do mkdir -p $i/trim_data;done
ls -d SRR107526*|sort|while read i;do mv $i/*.fq $i/trim_data;done 
#Mapping of Reads with Genome
#3.sh = hisat2  -x Tea.index -1 /home/scbb/abhijit/$1/trim_data/*_*_forward_paired.fq -2 /home/scbb/abhijit/$1/trim_data/*_*_reverse_paired.fq -S $1/out.sam -p 80
ls -d SRR107526* |while read i;do sh 3.sh ;done
#Quantification of Genes
ls -d SRR107526*|sort|while read i;do htseq-count  $i/*.sam  GCF_004153795.1_AHAU_CSS_1_genomic.gtf > $i/$i\_counts.txt;done
