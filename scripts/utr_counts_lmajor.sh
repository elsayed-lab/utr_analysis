#!/bin/env zsh

sl_subdir="minlength4_mindiff-2"
polya_subdir="minlength4_mindiff-2"

# spliced leader
for stage in "metacyclic" "procyclic"; do
    # top-level directory
    base_dir=$SCRATCH/utr_analysis/lmajor-${stage}

    echo "===================================="
    echo " Processing $stage (Spliced leader)"
    echo "===================================="

    # iterate over samples
    for x in ${base_dir}/common/*; do
        sample_id=$(basename $x)

        # directories
        sldir="${base_dir}/spliced_leader/${sl_subdir}/${sample_id}/results"
        rsldir="${base_dir}/reverse_spliced_leader/${sl_subdir}/${sample_id}/results"

        # sub-totals
        slr1=$(tail -n+5 ${sldir}/sl_coordinates_R1.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        slr2=$(tail -n+5 ${sldir}/sl_coordinates_R2.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        rslr1=$(tail -n+5 ${rsldir}/rsl_coordinates_R1.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        rslr2=$(tail -n+5 ${rsldir}/rsl_coordinates_R2.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)

        # total number of reads
        ((total=$slr1 + $slr2 + $rslr1 + $rslr2))

        printf "%s: %d\n" $sample_id $total
    done
done

# poly(a) sites
for stage in "metacyclic" "procyclic"; do
    # top-level directory
    base_dir=$SCRATCH/utr_analysis/lmajor-${stage}

    echo "===================================="
    echo " Processing $stage (PolyA)"
    echo "===================================="

    # iterate over samples
    for x in ${base_dir}/common/*; do
        sample_id=$(basename $x)

        # directories
        polyadir="${base_dir}/poly-a/${polya_subdir}/${sample_id}/results"
        polytdir="${base_dir}/poly-t/${polya_subdir}/${sample_id}/results"

        # sub-totals
        polyar1=$(tail -n+5 ${polyadir}/polya_coordinates_R1.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        polyar2=$(tail -n+5 ${polyadir}/polya_coordinates_R2.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        polytr1=$(tail -n+5 ${polytdir}/polyt_coordinates_R1.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)
        polytr2=$(tail -n+5 ${polytdir}/polyt_coordinates_R2.gff |\
                grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" | csvstat -c 6 -t --sum)

        # total number of reads
        ((total=$polyar1 + $polyar2 + $polytr1 + $polytr2))

        printf "%s: %d\n" $sample_id $total
    done
done
