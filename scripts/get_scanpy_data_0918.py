#!/usr/bin/env python 
import scipy.io as sio
import pandas as pd

'''
used on ratchet..
module load local_python/3.6.5 
'''

#Load Data
# data = sio.loadmat('20181203_SCH_80kcells_novaseq.mat')
# data = sio.loadmat('20181205_SCH_nextseq_80kcells_novaseq_qc.mat')
data = sio.loadmat('20190213_TimCherry_kit_test.mat')

#Digital Expression Matrix
DGE = data['DGE']

#Genes
genes = pd.Series(data['genes']).str.strip(' ')
print(genes)
# genes.to_csv('20181205_SCH_nextseq_80kcells_novaseq_qc.genes.csv')
genes.to_csv('20190213_TimCherry_kit_test.genes.csv')

#Sample types
sample_type = pd.Series(data['sample_type']).str.strip(' ')
# sample_type.to_csv('20181205_SCH_nextseq_80kcells_novaseq_qc.sample_type.csv')
sample_type.to_csv('20190213_TimCherry_kit_test.sample_type.csv')
# print(sample_type)


# #Main cluster assignment
# cluster_assignment = pd.Series(data['cluster_assignment']).str.strip(' ')

# #Spinal cluster assignment
# spinal_cluster_assignment = pd.Series(data['spinal_cluster_assignment']).str.strip(' ')


# sio.mmwrite('20181205_SCH_nextseq_80kcells_novaseq_qc.mtx', DGE)
sio.mmwrite('20190213_TimCherry_kit_test.mtx', DGE)