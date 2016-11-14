#!/bin/bash

###Script to Process Count RMEAN values from virus extracted data from Gene Omnibus####
#######################################################################################

##Extracting Human specific EBOLA virus related data for further analysis#################
##########################################################################################

cat EBOLA_RELATED_VIRUSES.txt | egrep 'RMEAN|Homo sapiens' > EBOLA_Human_related_viruses.txt

egrep -w -f virus_list.txt EBOLA_Human_related_viruses.txt | cut -f2 | sort | uniq > Samples_titles_used.txt

###Creating Title Specific files for gene and RMEAN count data Processing################
#########################################################################################

cat Samples_titles_used.txt | awk '{print "\"" $0 "\"", "EBOLA_Human_related_viruses.txt >"}' | sed s'/ \"/\"/'g | awk '{print "egrep -w", $0}' | nl | awk '{print $0, "sample_Iden_"$1".mat"}' | cut -f2 | awk '{gsub(/\(/,"\\(");gsub(/\)/,"\\)");print}' > Extract_All_sample_info.txt

echo "\#\!\/bin\/bash" | cat - Extract_All_sample_info.txt | awk '{gsub(/\\/,"");print}' > Extract_All_sample_info.txt-1 && mv Extract_All_sample_info.txt-1 Extract_All_sample_info.txt

chmod 755 Extract_All_sample_info.txt

sh Extract_All_sample_info.txt  && rm  Extract_All_sample_info.txt

##Ordering Virus title files for data extraction######
######################################################

find ./ -name '*.mat' | xargs grep -c ^.* | tr ':' '\t' | awk '{gsub(/\.\/\//,"");print}' | sort -k2,2nr | cut -f1 | grep 'sample_Iden' > Order_sample_names.txt

cp Order_sample_names.txt Samples_used.txt

###Extracting all Ebola Gene Names with gene expression variation from micro array data########
###############################################################################################

##Signaling specific genes expressed in Zaire Ebola data##################### 

echo "MVP POLG HLA SRC LCK MERTK TYRO3 UFO ACK1 BMPR2 YES FYN ROCK2 PAK1 TBK1 TTBK1 IKKE IMB1 IKBA IL1 IL2 IL3 IL4 IL5 IL6 IL7 IL8 IL9 IL10 CD209 MYB" | tr ' ' '\n' > Names.txt

##Genes down regulated in Zaire Ebola virus data#### 

cat EBOLA_Human_related_viruses.txt | egrep 'Zaire ebolavirus' | awk '$10<=25' | sort -k4,4n | cut -f5 | sed s'/<:>/:/'g | tr ':' '\n' | sed '/Homo/d' | sort | uniq | head -n 150 > List-1.txt

##Genes up regulated in Zaire Ebola virus data#### 

cat EBOLA_Human_related_viruses.txt | egrep 'Zaire ebolavirus' | awk '$10>=75' | sort -k4,4nr | cut -f5 | sed s'/<:>/:/'g | tr ':' '\n' | sed '/Homo/d' | sort | uniq | head -n 150 > List-2.txt 

cat List-1.txt List-2.txt Names.txt | sort | uniq > Ebola_Gene_Names.txt

rm List-1.txt List-2.txt Names.txt

##Obtaining Total gene count sum for Zaire Ebola virus data##

cat Ebola_Gene_Names.txt | wc -l > protein_count.txt

cat protein_count.txt | awk '{print $1}' > File.txt && mv File.txt protein_count.txt

##Running Extraction script##

chmod 755 Geo_profile_script.sh

sh Geo_profile_script.sh 

##Normalized RMEAN Count values###

mv counts_ALL.txt-1 counts_ALL-one.txt

##Raw RMEAN counts#################

mv Group_ALL.txt-1 Group_ALL-one.txt

mv counts_ALL-one.txt counts_ALL.txt

mv Group_ALL-one.txt Group_ALL.txt

###Putting a header for gene names##

echo "Gene_Name" | cat - Ebola_Gene_Names.txt > File.txt && mv File.txt Ebola_Gene_Names.txt

###attaching genes to their relative count data###

##Normalized data####

paste Ebola_Gene_Names.txt counts_ALL.txt | awk '{print $0}' | tr ' ' '\t' > Table_genes.txt

##Raw count data########

paste Ebola_Gene_Names.txt Group_ALL.txt | awk '{print $0}' | tr ' ' '\t' > CountTable_genes.txt

###Ensuring the data datable has only gene related data###################
##########################################################################

##Normalized Table

cat Table_genes.txt | grep '^[A-Z]' > File.txt && mv File.txt Table_genes.txt

###Raw count Table####

cat CountTable_genes.txt | grep '^[A-Z]' > CountTable_genes.txt-1 && mv CountTable_genes.txt-1 CountTable_genes.txt

rm *counts* *Group_ALL*

###Creating Headers for Each set of Table values####

cp Table_genes.txt Modified_table-genes.txt

cat Table_genes.txt | sed '1d' > File.txt && mv File.txt Table_genes.txt

cat CountTable_genes.txt | sed '1d' > File.txt && mv File.txt CountTable_genes.txt

cat Modified_table-genes.txt | sed -n '1p' | tr '\t' '\n' | sed '1d' | nl | awk '{gsub(/Sample_/,"Sample\t");print}' | awk '{print $2"_" $1, $3}' | tr ' ' '\t' | cut -f1 | tr '\n' ' ' | awk '{gsub(/ /,"\t");print}' | awk '{print "echo", "\"" "Gene_Name\t" $0"\""}' | awk '{gsub(/\\t\"/,"\"");print}' | awk '{print $0, "| cat - Table_genes.txt > Table_genes.txt-1 && mv Table_genes.txt-1 Table_genes.txt"}' > Table_header.txt

echo "\#\!\/bin\/bash" | cat - Table_header.txt | awk '{gsub(/\\/,"");print}' > File.txt && mv File.txt Table_header.txt

chmod 755 Table_header.txt

###Adding Sample headers to data table######

sh Table_header.txt && rm Table_header.txt

cat Modified_table-genes.txt | sed -n '1p' | tr '\t' '\n' | sed '1d' | nl | awk '{gsub(/Sample_/,"Sample\t");print}' | awk '{print $2"_" $1, $3}' | tr ' ' '\t' | cut -f1 | tr '\n' ' ' | awk '{gsub(/ /,"\t");print}' | awk '{print "echo", "\"" "Gene_Name\t" $0"\""}' | awk '{gsub(/\\t\"/,"\"");print}' | awk '{print $0, "| cat - CountTable_genes.txt > CountTable_genes.txt-1 && mv CountTable_genes.txt-1 CountTable_genes.txt"}' > Table_header.txt

echo "\#\!\/bin\/bash" | cat - Table_header.txt | awk '{gsub(/\\/,"");print}' > File.txt && mv File.txt Table_header.txt

chmod 755 Table_header.txt

sh Table_header.txt && rm Table_header.txt

###Creating a specific data table with LCK gene and related signaling protein expressions###### 
###############################################################################################

cat Table_genes.txt | egrep -w 'Gene_Name|MVP|POLG|HLA|SRC|LCK|MERTK|TYRO3|UFO|ACK1|BMPR2|YES|FYN|ROCK2|PAK1|TBK1|TTBK1|IKKE|IMB1|IKBA|CD209|MYB|IL2|IL1A|IL1B|IL1R1|IL1R2|IL1RAP|IL1RAPL1|IL1RAPL2|IL1RL1|IL1RL2|IL1RN|IL2RA|IL2RB|IL2RG' > Top_4_Ebola_Related_genes.txt

###Ensuring the data datable has only gene related data###################
##########################################################################

cat CountTable_genes.txt | egrep '^[A-Z]' > Counts_Ebola_Related_genes.txt && mv Counts_Ebola_Related_genes.txt CountTable_genes.txt

mv Ebola_Gene_Names.txt Ebola_Protein_names.txt

echo "That was a great job"

echo "well done, You can now analyze your table \n For your desired expected outcome"

rm *.mat


