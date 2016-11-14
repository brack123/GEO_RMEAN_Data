#!/bin/bash

####Script For Extracting EBOLA RMEAN data per viral sample##################
#############################################################################

##Running the first 10 samples###

cat Order_sample_names.txt | sed -n '1,10p' > Group-S_sample.txt 

cat Order_sample_names.txt | sed '1,10d' > File.txt && mv File.txt Order_sample_names.txt  

###Designing the extraction script#########################################
###########################################################################

cat Group-S_sample.txt | nl | awk '{print "for i in `cat Ebola_Gene_Names.txt`; do egrep -w $i", $2, "| cut -f4 | uniq > $i-"$1";", "done"}' > Group-S_sample.txt-1 && mv Group-S_sample.txt-1 Group-S_sample.txt

cat Group-S_sample.txt | awk '{gsub(/-10/,"-0");print}' > Group-S_sample.txt-1 && mv Group-S_sample.txt-1 Group-S_sample.txt

echo "\#\!\/bin\/bash" | cat - Group-S_sample.txt | awk '{gsub(/\\/,"");print}' > File.txt && mv File.txt Group-S_sample.txt

chmod 755 Group-S_sample.txt 

sh Group-S_sample.txt && rm Group-S_sample.txt

###Processing the Gene and count data extraction script######################

for i in ./*-1*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-1.txt && cp counts-1.txt group_1.txt

####Normalizing the count data across samples############

cat protein_count.txt | awk '{print "cat counts-1.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_1.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-1*

for i in ./*-2*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-2.txt && cp counts-2.txt group_2.txt

cat protein_count.txt | awk '{print "cat counts-2.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_2.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-2*

for i in ./*-3*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-3.txt && cp counts-3.txt group_3.txt

cat protein_count.txt | awk '{print "cat counts-3.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_3.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-3*

for i in ./*-4*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-4.txt && cp counts-4.txt group_4.txt

cat protein_count.txt | awk '{print "cat counts-4.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_4.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-4*

for i in ./*-5*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-5.txt && cp counts-5.txt group_5.txt

cat protein_count.txt | awk '{print "cat counts-5.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_5.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-5*

for i in ./*-6*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-6.txt && cp counts-6.txt group_6.txt

cat protein_count.txt | awk '{print "cat counts-6.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_6.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-6*

for i in ./*-7*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-7.txt && cp counts-7.txt group_7.txt

cat protein_count.txt | awk '{print "cat counts-7.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_7.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-7*

for i in ./*-8*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-8.txt && cp counts-8.txt group_8.txt

cat protein_count.txt | awk '{print "cat counts-8.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_8.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh 

sh script.sh && rm *-8*

for i in ./*-9*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-9.txt && cp counts-9.txt group_9.txt

cat protein_count.txt | awk '{print "cat counts-9.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_9.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh

sh script.sh && rm *-9*

for i in ./*-0*; do echo "NA:0" | cat - $i | awk '{print "S" $0}'| tr '\n' '\t' | awk '{gsub(/:SNA/,"\nSNA");print}' | sed s'/:0:/:/'g | awk '{gsub(/SNA:0\t/,"");print}' | awk '{gsub(/\t/,"+0");print}' | sed s'/S//'g | awk '{print "0"$0}' | bc; done | cat > counts-0.txt && cp counts-0.txt group_0.txt

cat protein_count.txt | awk '{print "cat counts-0.txt |", "awk \"{print", "\\""$1,", "\\""\""$0"\\""\"""}\"", "| awk \"{print", "\\""$1/""\\""$2""}\"", "> counts_0.txt"}' > script.sh 

cat script.sh | awk '{print "#!/bin/bash:"$0}' | tr ':' '\n' > script.sh-1 && mv script.sh-1 script.sh

chmod 755 script.sh  

sh script.sh && rm *-0*

paste counts_1.txt counts_2.txt counts_3.txt counts_4.txt counts_5.txt counts_6.txt counts_7.txt counts_8.txt counts_9.txt counts_0.txt > counts_ALL.txt

paste group_1.txt group_2.txt group_3.txt group_4.txt group_5.txt group_6.txt group_7.txt group_8.txt group_9.txt group_0.txt > Group_ALL.txt

echo "BC_Sample_1\tBC_Sample_2\tBC_Sample_3\tBC_Sample_4\tBC_Sample_5\tBC_Sample_6\tBC_Sample_7\tBC_Sample_8\tBC_Sample_9\tBC_Sample_10" |  cat - counts_ALL.txt > counts_ALL.txt-1 

echo "BC_Sample_1\tBC_Sample_2\tBC_Sample_3\tBC_Sample_4\tBC_Sample_5\tBC_Sample_6\tBC_Sample_7\tBC_Sample_8\tBC_Sample_9\tBC_Sample_10" |  cat - Group_ALL.txt > Group_ALL.txt-1 

rm *counts_[0-9]* *group_[0-9]*





