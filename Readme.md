# vcf_annotator (DNAnexus Platform App)

Using one source vcf, annotate another vcf.

Example use:
```
dx run vcf_annotator -iraw_vcf=raw.vcf -isrc_vcf=src.vcf -ifields="ID,QUAL,+TAG" -o annotated.vcf.gz
```

Make sure the src_vcf is bgzipped using:
- `bgzip src_vcf`
- `bcftools view src_vcf -Oz`