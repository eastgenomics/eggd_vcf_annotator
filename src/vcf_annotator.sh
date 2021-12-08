#!/bin/bash

set -euxo pipefail

main() {

    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of src_vcf: '$src_vcf'"
    echo "Value of fields: '$fields'"
    if [ -z ${output_suffix+x} ]; then
        output_suffix="custom_annotated"
    fi

    echo "Value of output_suffix: '$output_suffix'"

    dx download "$raw_vcf"
    dx download "$src_vcf"

    # Compile bcftools
    tar xvjf bcftools-1.14.tar.bz2

    cd bcftools-1.14

    ./configure
    make
    make install

    cd ..

    # Bgzip the raw vcf
    bcftools view $raw_vcf_name -Oz > ${raw_vcf_name}.gz

    # Index both given vcfs
    bcftools index ${raw_vcf_name}.gz
    bcftools index ${src_vcf_name}

    # define file output name
    basename=$(echo $raw_vcf_name | cut -d"." -f1)
    annotated_vcf_file=${basename}_${output_suffix}.vcf
    annotated_vcf_bgzip_file=${annotated_vcf_file}.gz

    # add prefix to fields that need them
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

    # Annotate and bgzip output
    bcftools annotate -a $src_vcf_name -c $new_fields ${raw_vcf_name}.gz > $annotated_vcf_file
    bcftools view ${annotated_vcf_file} -Oz > $annotated_vcf_bgzip_file

    annotated_vcf=$(dx upload $annotated_vcf_bgzip_file --brief)

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
}
