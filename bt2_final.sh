#!/bin/bash

# Start full logging
mkdir -p logs
exec > >(tee logs/full_pipeline.log) 2>&1

# Create output directories
mkdir -p ./bt_output

# Loop through all .fastq.gz files
for file in ./*.fastq.gz; do
    base_name=$(basename "$file" .fastq.gz)
    output_file="./temp.fastq"
    echo ""
    echo ""
    echo "🔬 Processing: $base_name"

    # Trimming and filtering with fastp
    fastp --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
          --low_complexity_filter --length_required 70 --thread 16 \
          -i "$file" -o "$output_file"
          
    echo ""
    echo "fastp analysis of $file ended"
    echo ""
    
    # Taxonomic profiling / Alignment for each of 28 indices
    for i in {01..28} ; do 
        idx="/media/radmir/4A981D6C2CDB08E3/bacgen/bt2_${i}_index"
        out="./bt_output/${base_name}_${i}"

        # Align using Bowtie2 and convert to BAM
        echo ""
        echo "bt2 analysis for index $i"
        echo ""

        bowtie2 --sensitive --threads 32 -k 1 -x "$idx" -U "$output_file" | \
        samtools view -@ 32 -b -F 4 -q 30 -u - | \
        samtools sort -@ 32 -T "${out}_temp" -o "${out}.bam"


        # Index the BAM
        samtools index -@ 32 "${out}.bam"

        # Output idxstats
        samtools idxstats "${out}.bam" > "${out}.idxstats"

        # Generate per-read alignment coverage (%)
        samtools view -@ 32 -h -F 4 "${out}.bam" | awk '
        {
          if ($0 ~ /^@/) next;
          ref = $3;
          seq_len = length($10);
          cigar = $6;
          aligned = 0;

          while (match(cigar, /[0-9]+M/)) {
            val = substr(cigar, RSTART, RLENGTH);
            sub(/M/, "", val);
            aligned += val;
            cigar = substr(cigar, RSTART + RLENGTH);
          }

          if (seq_len > 0) {
            percent = 100 * aligned / seq_len;
            printf "%s\t%.4f\n", ref, percent;
          }
        }' > "${out}_alignment_coverage_raw.tsv"

        # Combine average percent with idxstats
        awk 'FNR==NR {sum[$1]+=$2; count[$1]++; next} 
             {
               avg = (count[$1] > 0) ? sum[$1]/count[$1] : 0;
               print $1, $2, $3, $4, avg;
             }' "${out}_alignment_coverage_raw.tsv" "${out}.idxstats" \
             > "${out}_alignment_coverage_summary.tsv"

                # Index the sorted BAM and output idxstats
        #        samtools index -@ 64 "${out}.sorted.bam"
        #        samtools idxstats "${out}.sorted.bam" > "${out}.idxstats"
    done



    # Clean up temporary fastp files
    echo ""
    echo "deleting fastp outputs, bam and bai files"
    rm -rf "$output_file" fastp.html fastp.json ./bt_output/*.bam ./bt_output/*.bam.bai ./bt_out/*.sam

    echo ""
    echo "✅ Done: $base_name"
    echo "----------------------------------------"
done

#s
seq_list=(10 12 13 14)# 15 16 17 18 19 1 20 21 22 23 24 25 26 2 32 33 37 38 39 3 40 41 42 43 44 45 46 47 48 49 4 50 51 52 53 54 55 56 57 58 59 5 61 62 63 64 7 9)

#c
seq_list=(s10 s12 s13 s14 s15 s16 s17 s18 s19 s1 s20 s21 s22 s23 s24 s25 s26 s2 s32 s33 s37 s38 s39 s3 s40 s41 s42 s43 s44 s45 s46 s47 s48 s49 s4 s50 s51 s52 s53 s54 s55 s56 s57 s58 s59 s5 s7 s9)

#fg
seq_list=(s10 s12 s13 s14 s15 s16 s17 s18 s19 s1 s20 s21 s22 s23 s24 s25 s26 s2 s32 s33 s37 s38 s39 s3 s40 s41 s42 s43 s44 s45 s46 s47 s48 s49 s4 s50 s51 s52 s53 s54 s55 s56 s57 s58 s59 s5 s7 s9)

for i in "${seq_list[@]}"; do
    echo "📦 Compressing: $i"
    tar -czf "tar.${i}.tar.gz" ./bt_output/${i}* || { echo "❌ Compression of $i failed."; exit 1; }

    # Verify archive integrity
    if tar -tzf "tar.${i}.tar.gz" > /dev/null 2>&1; then
        echo "🧹 Deleting original files starting with '$i'..."
        rm -v ./bt_output/${i}*
        echo "✅ Done."
    else
        echo "❌ Integrity check failed for tar.${i}.tar.gz. Files not deleted."
    fi
done