
#get protein seq

../datasets download gene 59272 70008 492331 712790 112313373 --filename ACE2.zip
unzip ACE2.zip -d ACE2
#ls ACE2/ncbi_dataset/data/protein.faa

../datasets download gene 28, 80908, 722252, 112320051 --filename ABO.zip
unzip ABO.zip -d ABO
#ls ABO/ncbi_dataset/data/protein.faa

../datasets download gene 7113, 50528, 494080, 715138, 112306012 --filename TMPRSS2.zip
unzip TMPRSS2.zip -d TMPRSS2
#ls TMPRSS2/ncbi_dataset/data/protein.faa


#extract regions of protein seq

echo -n "" > tempRegion

for i in $(grep '>'  ./ABO/ncbi_dataset/data/protein.faa | cut -f 1 -d ' '| sed 's/>//')
do
        echo $i
        #esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}'
        esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd sig_peptide   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "sig_peptide\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd mat_peptide   product INSDInterval_from INSDInterval_to sub_sequence | awk '{print "mat_peptide\t"$0}' >> tempRegion

done

paste tempRegion <(cat tempRegion| cut -f 3 | cut -f 1 -d '/') | awk  'BEGIN {FS="\t"}{print ">"$2"|"$1"|"$4"-"$5"|"$8"\n"$6}' > ABO_feature.fasta


echo -n "" > tempRegion

for i in $(grep '>'  ./ACE2/ncbi_dataset/data/protein.faa | cut -f 1 -d ' '| sed 's/>//')
do
        echo $i
        esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd sig_peptide   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "sig_peptide\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd mat_peptide   product INSDInterval_from INSDInterval_to sub_sequence | awk '{print "mat_peptide\t"$0}' >> tempRegion

done

paste tempRegion <(cat tempRegion| cut -f 3 | cut -f 1 -d '/') | awk  'BEGIN {FS=OFS="\t"}{print ">"$2"|"$1"|"$4"-"$5"|"$8"\n"$6}' > ACE2_feature.fasta


echo -n "" > tempRegion

for i in $(grep '>'  ./TMPRSS2/ncbi_dataset/data/protein.faa | cut -f 1 -d ' '| sed 's/>//')
do
        echo $i
        esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd sig_peptide   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "sig_peptide\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd mat_peptide   product INSDInterval_from INSDInterval_to sub_sequence | awk '{print "mat_peptide\t"$0}' >> tempRegion

done

paste tempRegion <(cat tempRegion| cut -f 3 | cut -f 1 -d '/') | awk  'BEGIN {FS=OFS="\t"}{print ">"$2"|"$1"|"$4"-"$5"|"$8"\n"$6}' > TMPRSS2_feature.fasta

#ls ABO_feature.fasta ACE2_feature.fasta TMPRSS2_feature.fasta
