#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA as sklearnPCA
import argparse as ap

#parse arguments
parser = ap.ArgumentParser()
parser.add_argument('-i', '--infile', help='A haplotypes.tsv file as created by populations in the Stacks package (Catchen et al 2013)', required=True, type=str)
parser.add_argument('-n', '--max_snps', help='Maximum number of SNPs per locus', required=True, type=int)
parser.add_argument('-m', '--max_missing_data', help='Maximum %% of missing loci to be included in the PCA', required=True, type=int)

args = parser.parse_args()

#set input and output (not input of fineRADpainter) file
input_file = open(args.infile, "r")
output_file = args.infile
output_file += '.lociFilt.txt'
output = open(output_file, 'w')

#variables
max_snps = args.max_snps

#set counters
triallelic = 0
invariants = 0
too_many_snps = 0
indels = 0

#set lists
num_alleles_list = []
num_snp_list = []
cat_id_list = []

#filter loci
for line in input_file.readlines():
	#replace missing data value from "-" to ""
	line = line.replace('-','')
	row = line.rstrip().split('\t')
	if row[0] == 'Catalog ID':
		header = row
		samples = header[2:len(header)]
		output.write(line)
	else:
		#separate genotypes into alleles
		row = [i.split('/') for i in row]
		#drop loci if more than 2 alleles are present (likely already filtered)
		if any(len(i) > 2 for i in row):
			triallelic += 1
		#drop invariant loci
		elif any(i == ["consensus"] for i in row):
			invariants += 1
			num_alleles_list.append(1)
			num_snp_list.append(0)
			cat_id_list.append(row[0])
		else:
			genodata = row[2:len(row)]
			all_alleles = [item for sublist in genodata for item in sublist]
			if any(len(ind) > max_snps for ind in all_alleles): #drop loci if they have more than max_snps allowed
				too_many_snps += 1
			#drop loci with "N" likely indels but we need to check
			elif any('N' in ind for ind in all_alleles):
				indels += 1

			else:
				alleles_set = set([allele for allele in all_alleles if allele != ''])
				num_alleles_list.append(len(alleles_set))
				num_snp_list.append(len(list(alleles_set)[0]))
				cat_id_list.append(row[0])
				output.write(line)
output.close()
input_file.close()

print 'Loci filtered sequentially as follows:'
print 'Triallelic =', triallelic, '; Invariants =', invariants, '; More than', max_snps, 'snps =', too_many_snps, '; Indels =', indels
print '---'
print 'Filtered loci for all samples saved to '+args.infile+'.filtered'
print '---'

#plot distribution of number of alleles across loci 
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot(111)
counts = np.bincount(num_alleles_list)
plt.bar(range(max(num_alleles_list)+1), counts, width=1, align='center', color='DarkKhaki')
plt.xticks(range(max(num_alleles_list)+1))
plt.xlim(-1, max(num_alleles_list)+1)
plt.ylabel('# LOCI')
plt.xlabel('# ALLELES')
plt.title(str(len(num_alleles_list))+' loci analyzed (including invariants)')
missing_plot_file = args.infile+'.alleles.pdf'
fig.savefig(missing_plot_file,bbox_inches='tight')
plt.close()
print 'Distribution of alleles per locus plotted in '+args.infile+'.alleles.pdf'
print '---'

#plot distribution of number of snps across loci
plt.ioff()
fig = plt.figure()
ax = fig.add_subplot(111)
counts = np.bincount(num_snp_list)
plt.bar(range(max(num_snp_list)+1), counts, width=1, align='center', color='DarkSeaGreen')
plt.xticks(range(max(num_snp_list)+1))
plt.xlim(-1, max(num_snp_list)+1)
plt.ylabel('# LOCI')
plt.xlabel('# SNPS')
plt.title(str(len(num_snp_list))+' loci analyzed (including invariants)')
missing_plot_file = args.infile+'.snps.pdf'
fig.savefig(missing_plot_file,bbox_inches='tight')
plt.close()
print 'Distribution of SNPs per locus plotted in '+args.infile+'.snps.pdf'
print '---'

#Read the output file to plot missing data per sample and PCA of missing data and prepare the fineRADpainter input files
data = pd.read_csv(output_file, sep='\t', header=0, usecols=samples)

fineRADpainter_input = args.infile+'.fineRADpainter.lociFilt.txt'
data.to_csv(fineRADpainter_input, header=True, sep='\t', index=False)
print 'FineRADpainter input file with filtered loci saved to '+args.infile+'.fineRADpainter.lociFilt.txt'
print '---'

num_row = data.shape[0]
missing_data = data.isnull().sum().apply(lambda x: round(x*1.0/num_row*100,1))


plt.ioff()
fig = plt.figure()
ax = fig.add_subplot(111)
missing_data.plot(kind='bar', ylim=(0,100), colormap='summer', ax = ax)
plt.ylabel('Missing data (%)')
plt.title(str(num_row)+' loci analyzed (excluding invariants)')
missing_plot_file = args.infile+'.missingdata.pdf'
fig.savefig(missing_plot_file,bbox_inches='tight')
plt.close()
print 'Missing data per sample plotted in '+args.infile+'.missingdata.pdf'
print '---'

fact = 100/args.max_missing_data

thr = data.shape[0]-data.shape[0]/fact
data_20_01 = data.dropna(thresh=thr, axis=1).notnull().astype('int')
sklearn_pca20 = sklearnPCA(n_components=2)
sklearn_pca20.fit(data_20_01)

filtered_missingX_file = args.infile+'.fineRADpainter.lociFilt.samples'+str(args.max_missing_data)+'%missFilt.txt'
data.dropna(thresh=thr, axis=1).to_csv(filtered_missingX_file, header=True, sep='\t', index=False)
print 'FineRADpainter input file with filtered loci and filtered samples with max missing loci =',str(args.max_missing_data),'% saved to '+args.infile+'.fineRADpainter.filt'+str(args.max_missing_data)+'%Missing'
print '---'

plt.ioff()
fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(sklearn_pca20.components_[0], sklearn_pca20.components_[1], 'ro' )

labels = list(data_20_01.columns.values)

for label, x, y in zip(labels, sklearn_pca20.components_[0] ,sklearn_pca20.components_[1] ):
    plt.annotate(label, xy=(x, y), xytext=(-20, 20), textcoords='offset points', ha='right', va='bottom', arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.title('PCA of missing data')
missing_plot_file = args.infile+'.missingdataPCA.pdf'
fig.savefig(missing_plot_file,bbox_inches='tight')
plt.close()
print 'PCA of missing data plotted in '+args.infile+'.missingdataPCA.pdf'
print '---'





