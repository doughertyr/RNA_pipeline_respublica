import sys
import numpy as np
import pandas as pd

out_file = sys.argv[1].split('/')[-1].rstrip('HaplotypeCaller_output.vcf') + 'AnnotateSNVs_output.txt'
print(out_file)

with open("{}".format(sys.argv[1], "r")) as f:
	vcf = f.readlines()[53:]
f.close()

chrom_vcf = []
snp = []
dbSNP = []
for l in vcf:
	l_split = l.split()
	chrom_vcf.append(l_split[0])
	snp.append(int(l_split[1]))
	if l_split[2] == '.':
		dbSNP.append(0)
	else:
		dbSNP.append(1)

vcf_data = np.column_stack((snp,chrom_vcf,dbSNP))
vcf_df = pd.DataFrame(vcf_data,columns=['snp','chrom','dbSNP'])
vcf_df.snp = vcf_df.snp.astype(int)


with open("ucsc_hg19_uniq.txt", "r") as f:
	ucsc = f.readlines()
f.close()

chrom_ucsc = []
gene = []
strand = []
tx_srt = []
tx_end = []
for l in ucsc:
	l_split = l.split()
	chrom_ucsc.append(l_split[0])
	gene.append(l_split[1])
	strand.append(l_split[2])
	tx_srt.append(l_split[3])
	tx_end.append(l_split[4])

data = np.column_stack((chrom_ucsc,gene,strand,tx_srt,tx_end))
ucsc_df = pd.DataFrame(data,columns=['chrom','gene','strand','tx_srt','tx_end'])

ucsc_df.tx_srt = ucsc_df.tx_srt.astype(int)
ucsc_df.tx_end = ucsc_df.tx_end.astype(int)

ucsc_df['length'] = abs(ucsc_df['tx_end'] - ucsc_df['tx_srt'])

CHROMOSOME = []
POS_0 = []
POS_1 = []
GENE = []
DBSNP = []
STRAND = []

chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
for i in chromosomes:
#for index, row in vcf_df[vcf_df['chrom'] == '1'].iterrows():
	for index, row in vcf_df[vcf_df['chrom'] == i].iterrows():
		try:
			G,S = ucsc_df.loc[ucsc_df[ucsc_df['chrom'] == 'chr{}'.format(i)][ucsc_df['tx_srt'] < row['snp']][ucsc_df['tx_end'] > row['snp']]['length'].idxmax()][['gene','strand']]
			#G,S = ucsc_df.loc[ucsc_df[ucsc_df['chrom'] == 'chr1'][ucsc_df['tx_srt'] < row['snp']][ucsc_df['tx_end'] > row['snp']]['length'].idxmax()][['gene','strand']]
			if G[:2] == "NR":
				GENE.append('Inte')
				STRAND.append('#')
			else:
				GENE.append(G)
				STRAND.append(S)

		except ValueError:
			GENE.append('Inte')
			STRAND.append('#')

		CHROMOSOME.append('chr{}'.format(i))
		#CHROMOSOME.append('chr1')
		POS_0.append(row['snp'])
		POS_1.append(row['snp']+1)
		DBSNP.append(row['dbSNP'])

col = ['chr','pos-0','pos-1','gene','dbSNP','strand']

with open(out_file,'w') as w:
	for i in range(len(GENE)):
		w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(CHROMOSOME[i],POS_0[i],POS_1[i],GENE[i],DBSNP[i],STRAND[i]))
w.close()
