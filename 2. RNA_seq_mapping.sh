#!/bin/bash #shebang:it tells the shell what program to interpret the script with, when executed
## Mapping reads to the reference genenome
### for loop from /star_folder/leaf_trimmed


for i in '/star_folder/leaf_trimmed/A632.P1_1.fastq.gz /star_folder/leaf_trimmed/A632.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/A632.P2_1.fastq.gz /star_folder/leaf_trimmed/A632.P2_2.fastq.gz' \
'/star_folder/leaf_trimmed/B73.1_1.fastq.gz /star_folder/leaf_trimmed/B73.1_2.fastq.gz' \
'/star_folder/leaf_trimmed/B73.2_1.fastq.gz /star_folder/leaf_trimmed/B73.2_2.fastq.gz' \
'/star_folder/leaf_trimmed/B73.3_1.fastq.gz /star_folder/leaf_trimmed/B73.3_2.fastq.gz' \
'/star_folder/leaf_trimmed/CML91.P1_1.fastq.gz /star_folder/leaf_trimmed/CML91.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/CML91.P2_1.fastq.gz /star_folder/leaf_trimmed/CML91.P2_2.fastq.gz' \
'/star_folder/leaf_trimmed/F7.P1_1.fastq.gz /star_folder/leaf_trimmed/F7.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/F7.P2_1.fastq.gz /star_folder/leaf_trimmed/F7.P2_2.fastq.gz' \
'/star_folder/leaf_trimmed/H99.1_1.fastq.gz /star_folder/leaf_trimmed/H99.1_2.fastq.gz' \
'/star_folder/leaf_trimmed/H99.2_1.fastq.gz /star_folder/leaf_trimmed/H99.2_2.fastq.gz' \
'/star_folder/leaf_trimmed/H99.3_1.fastq.gz /star_folder/leaf_trimmed/H99.3_2.fastq.gz' \
'/star_folder/leaf_trimmed/HP301.P1_1.fastq.gz /star_folder/leaf_trimmed/HP301.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/HP301.P2_1.fastq.gz /star_folder/leaf_trimmed/HP301.P2_2.fastq.gz' \
'/star_folder/leaf_trimmed/Mo17.P1_1.fastq.gz /star_folder/leaf_trimmed/Mo17.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/Mo17.P2_1.fastq.gz /star_folder/leaf_trimmed/Mo17.P2_2.fastq.gz' \
'/star_folder/leaf_trimmed/W153.P1_1.fastq.gz /star_folder/leaf_trimmed/W153.P1_2.fastq.gz' \
'/star_folder/leaf_trimmed/W153.P2_1.fastq.gz /star_folder/leaf_trimmed/W153.P2_2.fastq.gz' \;
do
$RUN STAR --runThreadN 18 \
      --sjdbGTFfile /glab/references/zea_mays/annotation/Zea_mays.B73_RefGen_v4.45.gtf \
      --genomeDir /genomes/B73/B73_star_index \
      --readFilesIn $i \
      --readFilesCommand zcat \
      --outFileNamePrefix /star_folder/B73_mapping/alignments.leaf_STAR/$(basename -z -s .fastq.gz $i) \
      --outSAMtype BAM SortedByCoordinate \
      --outMultimapperOrder Random \
      --outBAMsortingThreadN 2 \
      --outWigType bedGraph \
      --outWigStrand Stranded \
      --quantMode GeneCounts \
      --twopassMode Basic \
      --outReadsUnmapped Fastx \
      --outWigNorm RPM \
      --limitOutSJcollapsed 2000000
done
