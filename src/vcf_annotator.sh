#!/bin/bash

set -euxo pipefail

main() {

    echo "Value of raw_vcf: '$raw_vcf'"
    echo "Value of src_vcf: '$src_vcf'"
    echo "Value of reference_genome: '$reference_genome'"
    echo "Value of fields: '$fields'"

    if [ -z ${output_suffix+x} ]; then
        output_suffix="custom_annotated"
    fi

    echo "Value of output_suffix: '$output_suffix'"

    dx download "$raw_vcf"
    dx download "$src_vcf"
    dx download "$reference_genome"

    # get nb of cpus
    nb_cpus=$(grep -c ^processor /proc/cpuinfo)

    # Bgzip the raw vcf
    bcftools view --threads $nb_cpus $raw_vcf_name -Oz > ${raw_vcf_name}.gz

    # Index both given vcfs
    bcftools index --threads $nb_cpus ${raw_vcf_name}.gz
    bcftools index --threads $nb_cpus ${src_vcf_name}

    # define file output name
    basename=$(echo $raw_vcf_name | cut -d"." -f1)
    annotated_vcf_file=${basename}_${output_suffix}.vcf

    # add EGGD_ prefix to fields unless explicitly renamed otherwise
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

    # decompose/normalise/left align raw vcf
    bcftools norm --threads $nb_cpus -f $reference_genome_name -m -any $raw_vcf_name > decom_norm_raw.vcf

    # Annotate and bgzip output
    bcftools annotate --threads $nb_cpus -a $src_vcf_name -c $new_fields decom_norm_raw.vcf > decom_norm_raw_annoated.vcf

    # join multi allelic back
    #mv decom_norm_raw_annoated.vcf $annotated_vcf_file
    bcftools norm --threads $nb_cpus -f $reference_genome_name -m +any decom_norm_raw_annoated.vcf > $annotated_vcf_file

    annotated_vcf=$(dx upload $annotated_vcf_file --brief)

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
}
