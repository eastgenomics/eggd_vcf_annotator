<!-- dx-header -->
# vcf_annotator (DNAnexus Platform App)

## What does this app do?
Using one source vcf, annotate another vcf.

## What are typical use cases for this app?
This app may be executed as a standalone app.

## What data are required for this app to run?
This app requires 2 VCFs and a string to run.  
The first vcf is the one you want to annotate, the second is the VCF containing the annotation. The fields input is comma-separated string to specify which columns to annotate with.  
Be sure that the src_vcf is bgzipped using:
- `bgzip src_vcf`
- `bcftools view src_vcf -Oz`

## What does this app output?
This app outputs a bgzipped VCF file.

Example cmd line:
```
dx run vcf_annotator -iraw_vcf=raw.vcf -isrc_vcf=src.vcf -ifields="ID,QUAL,+TAG" -o annotated.vcf.gz
```


#### This app was made by EMEE GLH
