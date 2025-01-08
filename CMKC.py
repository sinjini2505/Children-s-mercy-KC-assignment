import pysam
import pandas as pd

# read sample vcf
sample_file = "NA12878.chr21.slice.vcf"
sample_vcf = pysam.VariantFile(sample_file)
sample_headers = str(sample_vcf.header).splitlines()[-1][1:].split('\t')
sample_contents = []
for record in sample_vcf.fetch():
    sample_contents.append(str(record).split('\t'))
df_sample = pd.DataFrame(sample_contents, columns = sample_headers)

# read gnomad vcf
gnomad_file = "gnomad.chr21.slice.vcf"
gnomad_vcf = pysam.VariantFile(gnomad_file)
gnomad_headers = str(gnomad_vcf.header).splitlines()[-1][1:].split('\t')
gnomad_contents = []
for record in gnomad_vcf.fetch():
    gnomad_contents.append(str(record).split('\t'))
df_gnomad = pd.DataFrame(gnomad_contents, columns = gnomad_headers)

# filter gnomad vcf for population with allele frequency < `0.01`
df_gnomad['gnomad_AF'] = pd.to_numeric(df_gnomad.INFO.str.split(';').str[2].str.split('=').str[1])
df_gnomad = df_gnomad[df_gnomad['gnomad_AF']<0.01]

# variants in sample with gnomAD population allele frequency of less than `0.01`
df_merge = pd.merge(df_sample, df_gnomad, how='left', on=['CHROM', 'POS', 'REF', 'ALT'], suffixes=('_sample', '_gnomad'), indicator=True)
df_merge = df_merge[df_merge['_merge']=='both']

outfile = "output.xlsx"
df_merge.to_excel(outfile, index=False)