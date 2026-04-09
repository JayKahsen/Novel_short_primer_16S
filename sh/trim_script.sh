#!/bin/sh

clear

source /home/rushgmcf/miniconda3/bin/activate qiime2-2023.5

cd /home/rushgmcf/Jay/vplE2

# Section: Script Initialization
echo "Script 1"

# Section: Import Data
echo "Input metadata.tsv and manifest.tsv and .fastq"
qiime metadata tabulate --m-input-file metadata.tsv --o-visualization metadata.qzv
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./manifest.tsv \
  --output-path ./demux_seqs_original.qza

# Section: Original Data Visualization
echo "Visualize original data"
qiime demux summarize --i-data ./demux_seqs_original.qza --o-visualization ./demux_seqs_original.qzv

# Section: Trimming
echo "Trim"
qiime cutadapt trim-single \
  --i-demultiplexed-sequences demux_seqs_original.qza \
  --p-cores 8 \
  --p-minimum-length 240 \
  --o-trimmed-sequences trimmed_demux_seqs.qza \
  --p-adapter GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC \
  --p-discard-untrimmed
qiime demux summarize --i-data ./trimmed_demux_seqs.qza --o-visualization ./trimmed_demux_seqs.qzv

# Section: Denoising
# qiime dada2 denoise-paired more accurate?

echo "Denoise"
qiime dada2 denoise-single \
   --p-n-threads 10 \
   --i-demultiplexed-seqs ./trimmed_demux_seqs.qza \
   --p-trunc-len 0 \
   --output-dir DADA2_denoising_output \
   --verbose \
&> DADA2_denoising.log

# Section: Denoising Statistics and Representative Sequences
echo "Denoising stats"
qiime metadata tabulate --m-input-file DADA2_denoising_output/denoising_stats.qza --o-visualization denoising_stats.qzv

echo "Representative sequences"
qiime feature-table tabulate-seqs --i-data DADA2_denoising_output/representative_sequences.qza --o-visualization rep_seqs.qzv

# Section: Feature Table
echo "Feature table"
qiime feature-table summarize \
--i-table DADA2_denoising_output/table.qza \
--o-visualization feature_table.qzv \
--m-sample-metadata-file metadata.tsv

# Section: Classification of Representative Sequences
echo "Classify rep seqs"
qiime feature-classifier classify-sklearn \
  --i-classifier /home/rushgmcf/References/qiime2_2023_5/silva-138-99-nb-classifier.qza \
  --i-reads ./DADA2_denoising_output/representative_sequences.qza \
  --o-classification classified_rep_seqs.qza \
  --verbose \
  &> DADA2_classifier.log

# Section: Tabulation of Classified Sequences
echo "Tabulate the features, their taxonomy and the confidence of taxonomy assignment"
qiime metadata tabulate \
 --m-input-file classified_rep_seqs.qza \
 --o-visualization classified_rep_seqs.qzv

# Section: Export Feature Table to TSV
echo "Feature table to BIOM to TSV"
qiime tools export --input-path DADA2_denoising_output/table.qza --output-path DADA2_denoising_output
biom convert -i DADA2_denoising_output/feature-table.biom -o DADA2_denoising_output/collapsed_table.tsv --to-tsv

# Section: Script 2
echo "Script 2"
# Script 2

# Section: NO DECONTAMINATION STEPS
###### NO DECONTAMINATION STEPS ############

# Section: Filtering out Mitochondria and Chloroplast
echo "Filter out mitochiondria and chloroplast"
# Filter out mitochiondria and chloroplast
qiime taxa filter-table \
  --i-table DADA2_denoising_output/table.qza \
  --i-taxonomy classified_rep_seqs.qza \
  --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota \
  --o-filtered-table DADA2_denoising_output/table_filtered_v2.qza

# Section: Feature Table Summarization
echo "Feature table"
# Feature table
qiime feature-table summarize \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--o-visualization feature_table_filtered_v2.qzv \
--m-sample-metadata-file metadata.tsv

# Section: Filtering Sequences
echo "Filter sequences for those features that contained mito and chloro"
# Filter sequences for those features that contained mito and chloro
qiime taxa filter-seqs \
  --i-sequences DADA2_denoising_output/representative_sequences.qza \
  --i-taxonomy classified_rep_seqs.qza \
  --p-exclude Mitochondria,Chloroplast,Unassigned,Eukaryota \
  --o-filtered-sequences DADA2_denoising_output/representative_sequences_filtered_v2.qza

# Section: Representative Sequences Tabulation
echo "Representative sequences"
# Representative sequences
qiime feature-table tabulate-seqs \
--i-data DADA2_denoising_output/representative_sequences_filtered_v2.qza \
--o-visualization rep_seqs_filtered_v2.qzv

# Section: Classify Filtered Representative Sequences
echo "Classify rep seqs"
# # Classify rep seqs
qiime feature-classifier classify-sklearn \
  --i-classifier /home/rushgmcf/References/qiime2_2023_5/silva-138-99-nb-classifier.qza \
  --i-reads ./DADA2_denoising_output/representative_sequences_filtered_v2.qza \
  --o-classification classified_rep_seqs_filtered_v2.qza \
  --verbose \
  &> DADA2_classifier.log

# Section: Tabulation of Classified Filtered Sequences
echo "Tabulate the features, their taxonomy and the confidence of taxonomy assignment"
qiime metadata tabulate \
 --m-input-file classified_rep_seqs_filtered_v2.qza \
 --o-visualization classified_rep_seqs_filtered_v2.qzv

# Section: Export Filtered Feature Table to TSV
echo "FILTERED Feature table to BIOM to TSV"
qiime tools export --input-path DADA2_denoising_output/table_filtered_v2.qza --output-path ./
biom convert -i feature-table.biom -o DADA2_denoising_output/collapsed_table_filtered_v2.tsv --to-tsv

# Section: Taxonomic Analysis
echo "Taxonomic analysis"
qiime taxa barplot \
 --i-table DADA2_denoising_output/table_filtered_v2.qza \
 --i-taxonomy classified_rep_seqs_filtered_v2.qza \
 --m-metadata-file ./metadata.tsv \
 --o-visualization sample_barplots.qzv

# Section: Create Directory
mkdir alpha_rarefaction

#####################################################################################################
# rarefication 750
echo "750"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 750 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_750.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 750 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_750 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_750_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_750

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 750 --o-rarefied-table Rarefied_FeatureTable_750.qza
qiime tools export --input-path Rarefied_FeatureTable_750.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_750.tsv --to-tsv
#####################################################################################################

# Copy tables
cp DADA2_denoising_output/collapsed_table_filtered_v2.tsv FeatureTable_PostRemoval.tsv
cp DADA2_denoising_output/collapsed_table.tsv FeatureTable_NoRemoval.tsv
#####################################################################################################
#Taxa Tables
#####################################################################################################
echo "Pylum"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 2 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Phylum.tsv
rm -r taxtable

echo "Class"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 3 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Class.tsv
rm -r taxtable

echo "Order"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 4 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Order.tsv
rm -r taxtable

echo "Family"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 5 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Family.tsv
rm -r taxtable

echo "Genus"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 6 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Genus.tsv
rm -r taxtable

echo "Species"
qiime taxa collapse --i-table DADA2_denoising_output/table_filtered_v2.qza --i-taxonomy classified_rep_seqs_filtered_v2.qza --p-level 7 --output-dir taxtable
qiime tools export --input-path taxtable/collapsed_table.qza --output-path taxtable
biom convert -i taxtable/feature-table.biom -o taxtable/collapsed_table.tsv --to-tsv
cp taxtable/collapsed_table.tsv Table_Species.tsv
rm -r taxtable
#####################################################################################################
#####################################################################################################
# rarefication 1000
echo "1000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 1000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_1000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 1000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_1000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_1000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_1000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 1000 --o-rarefied-table Rarefied_FeatureTable_1000.qza
qiime tools export --input-path Rarefied_FeatureTable_1000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_1000.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 3000
echo "3000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 3000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_3000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 3000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_3000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_3000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_3000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 3000 --o-rarefied-table Rarefied_FeatureTable_3000.qza
qiime tools export --input-path Rarefied_FeatureTable_3000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_3000.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 5000
echo "5000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 5000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_5000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 5000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_5000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_5000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_5000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 5000 --o-rarefied-table Rarefied_FeatureTable_5000.qza
qiime tools export --input-path Rarefied_FeatureTable_5000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_5000.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 7500
echo "7500"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 7500 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_7500.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 7500 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_7500 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_7500_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_7500

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 7500 --o-rarefied-table Rarefied_FeatureTable_7500.qza
qiime tools export --input-path Rarefied_FeatureTable_7500.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_7500.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 10000
echo "10000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 10000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_10000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 10000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_10000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_10000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_10000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 10000 --o-rarefied-table Rarefied_FeatureTable_10000.qza
qiime tools export --input-path Rarefied_FeatureTable_10000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_10000.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 50000
echo "50000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 50000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_50000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 50000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_50000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_50000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_50000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 50000 --o-rarefied-table Rarefied_FeatureTable_50000.qza
qiime tools export --input-path Rarefied_FeatureTable_50000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_50000.tsv --to-tsv

#####################################################################################################
#####################################################################################################
# rarefication 100000
echo "100000"

# 16S rRNA dataset Alpha Diversity
qiime diversity alpha-rarefaction \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-max-depth 100000 \
--m-metadata-file ./metadata.tsv \
--p-metrics shannon --p-metrics simpson --p-metrics observed_features \
--o-visualization alpha_rarefaction/rarefaction_100000.qzv

qiime diversity core-metrics \
--i-table DADA2_denoising_output/table_filtered_v2.qza \
--p-sampling-depth 100000 \
--m-metadata-file ./metadata.tsv \
--output-dir core_metrics_100000 \
--p-n-jobs 8 \
--verbose \
&> core_metrics_100000_samples.log

# Run plugin for each alpha diversity result
cd core_metrics_100000

for result in *vector.qza; do
    outname=${result/_vector.qza/_group_significance.qzv}
    qiime diversity alpha-group-significance \
    --i-alpha-diversity $result \
    --m-metadata-file ../metadata.tsv \
    --o-visualization $outname
done

cd ..

# Feature table rarefaction
qiime feature-table rarefy --i-table DADA2_denoising_output/table_filtered_v2.qza --p-sampling-depth 100000 --o-rarefied-table Rarefied_FeatureTable_10000.qza
qiime tools export --input-path Rarefied_FeatureTable_100000.qza --output-path $PWD
biom convert -i feature-table.biom -o Rarefied_FeatureTable_100000.tsv --to-tsv

#####################################################################################################
