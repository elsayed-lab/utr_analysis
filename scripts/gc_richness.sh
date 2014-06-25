#!/bin/env bash
# GC Richness bedgraph generator
# 2014/04/27
#
# Source: http://wiki.bits.vib.be/index.php/Create_a_GC_content_track
#
# Usage example:
#
# gc_richness.sh /path/to/genome.fasta 10
genome=$1
width=$2

# generated files
genome_sizes="$(basename $genome).sizes"
genome_bed="$(basename $genome)_${width}bps.bed"
genome_freq="$(basename $genome)_nuc.txt"
output="${genome%.fasta}_GC_richness_${width}bps.igv"

# Generate list of chromosome lengths
samtools faidx $1
cut -f 1,2 ${genome}.fai > $genome_sizes

# Create bed file using specified window size
echo "Generating intermediate bed file"
bedtools makewindows -g $genome_sizes \
                     -w ${width} \
                     > $genome_bed

# Compute base frequencies
echo "Generating ${output}"
bedtools nuc -fi $genome -bed $genome_bed > $genome_freq

# Generate IGV formatted output file
echo "#track name=GC_richness color=83,124,166" > $output
gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"}
{
    if (FNR>1) {
        print $1,$2,$3,"GCpc_"w"bps",$5
    }
}' $genome_freq >> $output

# clean-up
rm $genome_sizes
rm $genome_bed
#rm $genome_freq

echo "Done!"

