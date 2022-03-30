import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

filename = '/home/noort/data/ref_files/tombo.DNA.txt'

df = pd.read_csv(filename,sep='\t', skiprows=1, header=None)
df.columns = ['kmer', 'mean', 'sd']
df.set_index('kmer', inplace = True)
seq = ''.join(np.random.choice(list('ACGT'), 1500))
df = df.sort_values('mean')
print(df)
print(seq)

plt.plot(df['mean'])
plt.show()
loc = np.arange(0, len(seq)-6)
sq = [df.loc[seq[i:i+6], 'mean'] for i in loc]
plt.plot(loc, sq)
plt.show()

tmp = 'GGTGTACTTCGTTCAGTTTACGTATTGCTAAGGTTAATAGGGAAACACGATACAGAATCCGGAGCAACTTTGTGTTTCCTCATCCACTCTCAATGCCTTCCTCAAGATCTGCTCAAAATTATTTTCATCTTTACCTGATGTTTCCTGTCACTCAGTTGGAGATGTTCCTCTGGCCATCAAGTCAGCCAGTTTAGTACCCCCTCTCCACCTTACTGGTACTCTCAAAGGGTGATTACCTACACCGTATGATTTTAATGATCTAAGCCCACAAGGCTCCTGAATCCAACTGACCTCCTCTCATCATCACTAAGGATTATTTCTTTGAGCTTGATGAATAGTCTTTCATTACATAATGCTGTCTATATAGAAATGTATGGATTTATATAAATTAAACTTGTGATTATTTCACTGATTTAGTCTATTGATTCAAACTAATCTTATGAGAGACATGAAAAAATAGTTCTTTTTGGTAGCTACTGACTGAGAGAGACAGTAGGGAATTTTCAGGGTCTGGAGGAACC'
print(len(tmp))