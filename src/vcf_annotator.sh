#!/bin/bash

set -euxo pipefail

main() {

    # Review input
    echo "Value of dest_vcf: '$dest_vcf'"
    echo "Value of src_vcf: '$src_vcf'"
    echo "Value of src_vcf_idx: '$src_vcf_idx'"
    echo "Value of ref_tar: '$ref_tar'"
    echo "Value of fields: '$fields'"

    # Append default output filename suffix if none provided
    # Avoids identical input/output vcf names
    if [ -z ${output_suffix+x} ]; then
        output_suffix="custom_annotated"
    fi

    echo "Value of output_suffix: '$output_suffix'"

    # Download inputs
    dx download "$dest_vcf"
    dx download "$src_vcf"
    dx download "$src_vcf_idx"
    dx download "$ref_tar"

    # Unpack reference genome tar
    tar xzf $ref_tar_name

    # Get nb of cpus to inform bcftools thread count
    nb_cpus=$(grep -c ^processor /proc/cpuinfo)

    # Add EGGD_ prefix to INFO fields unless explicitly renamed otherwise
    fields_array=()
    IFS=","
    for field in $fields; do
        if [[ $field != *":="* ]]; then
            field_to_add="EGGD_${field}:=${field}"
        else
            field_to_add=$field
        fi
        fields_array+=($field_to_add)
    done

    new_fields=$(IFS=","; echo "${fields_array[*]}")
    IFS=""

    # Decompress if input vcf is vcf.gz
    if [[ $dest_vcf_name == *.vcf.gz ]]; then
        zcat $dest_vcf_name > dest.vcf
    else
        mv $dest_vcf_name dest.vcf
    fi

    # Fix issue where vcfs have CSQ with no value, which causes bcftools norm
    # to fail. To fix this we change these to CSQ=.
    # CSQ is always the last field with input from nirvana2vcf so only
    # one sed pattern needed

    cat dest.vcf | \
    # CSQ is 1st INFO field
    sed s/"\\tCSQ;"/"\\tCSQ=.;"/g | \
    # CSQ is mid INFO field
    sed s/";CSQ;"/";CSQ=.;"/g | \
    # CSQ is last INFO field
    sed s/";CSQ\\t"/";CSQ=.\t"/g \
    > csq_fix.vcf
    bgzip csq_fix.vcf

    # Define file output name
    basename=$(echo $dest_vcf_name | cut -d"." -f1)
    annotated_vcf_file=${basename}_${output_suffix}.vcf.gz

    # Decompose/normalise/left align raw vcf, bgzip output
    bcftools norm --threads $nb_cpus -f genome.fa -m -any -Oz \
     csq_fix.vcf.gz > decom_norm_raw.vcf.gz
    bcftools index --threads $nb_cpus decom_norm_raw.vcf.gz

    # Annotate and bgzip output
    bcftools annotate --threads $nb_cpus -a $src_vcf_name -c $new_fields \
     decom_norm_raw.vcf.gz -Oz > decom_norm_raw_annotated.vcf.gz
    bcftools index --threads $nb_cpus decom_norm_raw_annotated.vcf.gz

    # Merge multi-allelic variants back into single records
    bcftools norm --threads $nb_cpus -f genome.fa -m +any \
     decom_norm_raw_annotated.vcf.gz -Oz > $annotated_vcf_file
    bcftools index --threads $nb_cpus --tbi $annotated_vcf_file

    # Note that the split and merge operation produces a minor change in the
    #  output vcf in addition to the requested annotation.
    #
    # PL for the called genotype is 0 in the input vcf, however for 
    #  multi-allelic calls this information is lost during the split/merge 
    #  operation since no single vcf record contains the called genotype 
    #  (i.e. 1/2) after the split, and so it cannot be restored in the merge
    #
    # Consequently PL of the called genotype becomes . in the output vcf

    annotated_vcf=$(dx upload $annotated_vcf_file --brief)
    annotated_vcf_index=$(dx upload ${annotated_vcf_file}.tbi --brief)

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
    dx-jobutil-add-output annotated_vcf_index "$annotated_vcf_index" --class=file
}
