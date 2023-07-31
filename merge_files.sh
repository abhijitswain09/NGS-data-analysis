#make a new directory
ls -d SRR107526*|sort|while read i; do mkdir -p $i/new ;done
#to create a column with all gene name
awk -F'\t' '{print $1}' SRR10752653/SRR10752653_counts.txt > merge1.txt
#extract the gene expression to a new folder
ls -d SRR107526*|sort|while read i;do awk '{print $2}' $i/$i\_counts.txt > $i/new/$i ;done
#merge the datas in a one file
cut -f1 merge1.txt | paste - SRR10752653/new/SRR10752653 SRR10752654/new/SRR10752654 SRR10752655/new/SRR10752655 SRR10752656/new/SRR10752656 SRR10752657/new/SRR10752657 SRR10752658/new/SRR10752658 SRR10752659/new/SRR10752659 SRR10752660/new/SRR10752660 SRR10752661/new/SRR10752661 SRR10752662/new/SRR10752662 SRR10752663/new/SRR10752663 SRR10752664/new/SRR10752664 SRR10752665/new/SRR10752665 SRR10752666/new/SRR10752666 SRR10752667/new/SRR10752667 > merge_all.txt
#print the header
awk  'BEGIN {print "Gene\tSRR10752653\tSRR10752654\tSRR10752655\tSRR10752656\tSRR10752657\tSRR10752658\tSRR10752659\tSRR10752660\tSRR10752661\tSRR10752662\tSRR10752663\tSRR10752664\tSRR10752665\tSRR10752666\tSRR10752667\tAdd"} {print $0"\t"$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16}' merge_all.txt > count_all_gene.txt
#delete the last 5 lines(30175-30179)
sed -i '30175,$d' count_all_gene.txt
#count the low expression dataset
awk -F'\t' '{if(!($17 <= 15)) print$0}' count_all_gene.txt |wc -l
 
 
