<!-- dx-header -->
# vcf_annotator (DNAnexus Platform App)

## What does this app do ?

Using one source vcf, annotate another vcf.

## What are typical use cases for this app ?

This app may be executed as a standalone app. One use case would be to add clinvar annotation to a patient vcf.

## What data are required for this app to run ?

This app requires 2 VCFs and a string to run.  
The first vcf is the one you want to annotate, the second is the VCF containing the annotation. The fields input follows the following snippet from the bcftools documentation (http://samtools.github.io/bcftools/bcftools.html#annotate):

> Comma-separated list of columns or tags to carry over from the annotation file (see also -a, --annotations). If the annotation file is not a VCF/BCF, list describes the columns of the annotation file and must include CHROM, POS (or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which should be ignored can be indicated by "-". If the annotation file is a VCF/BCF, only the edited columns/tags must be present and their order does not matter. The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags can be written both as "INFO/TAG" or simply "TAG", and FORMAT tags can be written as "FORMAT/TAG" or "FMT/TAG". The imported VCF annotations can be renamed as "DST_TAG:=SRC_TAG" or "FMT/DST_TAG:=FMT/SRC_TAG". To carry over all INFO annotations, use "INFO". To add all INFO annotations except "TAG", use "^INFO/TAG". By default, existing values are replaced. To add annotations without overwriting existing values (that is, to add missing tags or add values to existing tags with missing values), use "+TAG" instead of "TAG". To append to existing values (rather than replacing or leaving untouched), use "=TAG" (instead of "TAG" or "+TAG"). To replace only existing values without modifying missing annotations, use "-TAG". If the annotation file is not a VCF/BCF, all new annotations must be defined via -h, --header-lines. See also the -l, --merge-logic option.

Be sure that the src_vcf is bgzipped using:

- `bgzip src_vcf`
- `bcftools view src_vcf -Oz`

Example cmd line:

``` bash
dx run vcf_annotator -iraw_vcf=raw.vcf -isrc_vcf=src.vcf -ifields="ID,QUAL,+TAG" -o annotated.vcf.gz

dx run vcf_annotator -iraw_vcf=raw.vcf -isrc_vcf=src.vcf -ifields="TAG_RENAMED:=TAG" -o annotated.vcf.gz
```

## What does this app output?

This app outputs a bgzipped VCF file. The new fields in the vcf will be prefixed by "EGGD" unless specifically named otherwise in `fields` input.

Important note: Only `FILTER=PASS` variants from the src vcf will be annotated in the output vcf

### This app was made by EMEE GLH
