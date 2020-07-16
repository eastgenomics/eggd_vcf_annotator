#!/bin/bash

set -euxo pipefail

main() {

    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of src_vcf: '$src_vcf'"
    echo "Value of fields: '$fields'"

    dx download "$raw_vcf"
    dx download "$src_vcf"

    tar xvjf bcftools-1.10.2.tar.bz2
    
    cd bcftools-1.10.2

    ./configure
    make
    make install

    cd ..

    bcftools view $raw_vcf_name -Oz > ${raw_vcf_name}.gz

    bcftools index ${raw_vcf_name}.gz
    bcftools index ${src_vcf_name}

    basename=$(echo $raw_vcf_name | cut -d"." -f1)
    annotated_vcf_file=${basename}_annotated.vcf
    annotated_vcf_bgzip_file=${annotated_vcf_file}.gz

    bcftools annotate -a $src_vcf_name -c $fields ${raw_vcf_name}.gz > $annotated_vcf_file
    bcftools view ${annotated_vcf_file} -Oz > $annotated_vcf_bgzip_file

    annotated_vcf=$(dx upload $annotated_vcf_bgzip_file --brief)

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
}
