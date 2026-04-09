#!/bin/sh
cd $PWD
FILES=$PWD/*.gz
for f in $FILES
do
        fname=${f##*/}
        finalname4=${fname%_*}
        echo $finalname4
        mv "$fname" "$finalname4.fastq.gz"
done

#Pear Merging - merges and moves files merged reads, deletes other fastq files
mkdir PEAR_MERGEDFiles
echo -n "SampleName" >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo -n "Merged_Count" >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo -n "Raw_Count" >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo -n "Merging_Percentage" >> Pear_Summary.txt && echo -e ' \t ' >> Pear_Summary.txt
for f in *R1*.fastq.gz
do
        #echo $f
        name=$(basename ${f} _L001_R1.fastq.gz)
        #echo $name
        R2Name=${name}_L001_R2.fastq.gz
        #echo $R2Name
        Raw_Count=$(zcat $f | echo -n $((`wc -l`/4)))
        /home/rushgmcf/References/pear-0.9.11 -f $f -r $R2Name -o PEAR_MERGEDFiles/${name}.MERGED
        Merged_Count=$(cat PEAR_MERGEDFiles/${name}.MERGED.assembled.fastq | echo -n $((`wc -l`/4)))
        #echo $Merged_Count
        echo -n $name >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo -n $Merged_Count >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo -n $Raw_Count >> Pear_Summary.txt && echo -n -e ' \t ' >> Pear_Summary.txt && echo "$Merged_Count" "$Raw_Count" | awk '{print ($1)/$2*100}' >> Pear_Summary.txt
        rm $f
        rm $R2Name
done
mkdir demultiplexed_seqs
mv PEAR_MERGEDFiles/*.assembled.fastq demultiplexed_seqs/
rm -r PEAR_MERGEDFiles
rm *.fastq.gz
#Renaming
#MergedFiles=demultiplexed_seqs/*.fastq
cd demultiplexed_seqs/
for f in *.fastq
do
        fname=${f##*/}
        echo $fname
        trimname=${fname%_*}
	finalname31=${trimname//-/_}
       
        #trimname=${fname%.*}
        #finalname3=${trimname%.*}
        #finalname31=${finalname3%.*}
        #finalname33=${finalname31%_*}
        #finalname4=${finalname33%_*}
        echo $finalname31 >> ../Renaming_Output.txt
        mv "$fname" "$finalname31.fastq"
done
cd ..
tar -zcvf MergedFASTQ.tar.gz demultiplexed_seqs

