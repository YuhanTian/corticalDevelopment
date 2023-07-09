###############################################################################################################################################################
# run.scVelo_step1.R
###############################################################################################################################################################
# /home/yuchen/miniconda3/envs/R4.0/bin/R
library(Seurat)
library(Signac)
library(Matrix)
obj <- readRDS("/data3/yuhan/Project/Neuron/scMultiome/4_developmental_trajectories.newVersion.cleanedcell/1_URDTree/URDobj.rds")

# save metadata table:
obj$barcode <- colnames(obj)
obj$UMAP_1 <- obj@reductions$umap_combine@cell.embeddings[,1]
obj$UMAP_2 <- obj@reductions$umap_combine@cell.embeddings[,2]
write.csv(obj@meta.data, file='metadata.csv', quote=F, row.names=F)
# write expression counts matrix
counts_matrix <- GetAssayData(obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')
# write dimesnionality reduction matrix, in this example case harmony matrix
write.csv(obj@reductions$harmony@cell.embeddings, file='harmony.csv', quote=F, row.names=F)
# write gene names
write.table(data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',quote=F,row.names=F,col.names=F)


###############################################################################################################################################################
# run.scVelo_step2.py
#!/usr/bin/env python
# coding=utf-8

###############################################################################################################################################################
# conda env list
# conda create -n scVelo_update
conda activate scVelo_update
# pip install anndata
# pip install cellrank
ipython

import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("counts.mtx")
# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())
# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")
# load gene names:
with open("gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()
# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
# load dimensional reduction:
harmony = pd.read_csv("harmony.csv")
harmony.index = adata.obs.index
# set harmony and umap
adata.obsm['X_harmony'] = harmony.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
# plot a UMAP colored by Combine_CellType to test:
sc.pl.umap(adata, color=['Combine_CellType'], legend_loc='on data', legend_fontsize=10, legend_fontweight='normal', size=20, frameon=False, title='', save='UMAP.Combine_CellType.check1.pdf')
# save dataset as anndata format
# adata.write('AllCells_RAWadata.h5ad')
# reload dataset
# adata = sc.read_h5ad('AllCells_RAWadata.h5ad')


###############################################################################################################################################################
# run.scVelo_step3.py
#!/usr/bin/env python
# coding=utf-8

###############################################################################################################################################################
conda activate scVelo_update
ipython

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False, figsize=(7,7))
cr.settings.verbosity = 2

# reload dataset
# adata = sc.read_h5ad('AllCells_RAWadata.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('/data3/yuhan/Project/Neuron/scMultiome/scMultiome_library/4_GW20B1_FU_scMultiome_ATAC_GEX_B1/velocyto/4_GW20B1_FU_scMultiome_ATAC_GEX_B1_velocyto/gex_possorted_bam_6TT1Y.loom', cache=True)
ldata2 = scv.read('/data3/yuhan/Project/Neuron/scMultiome/scMultiome_library/8_GW15B1_301_scMultiome_ATAC_GEX_B1/velocyto/8_GW15B1_301_scMultiome_ATAC_GEX_B1_velocyto/gex_possorted_bam_PT2VL.loom', cache=True)
ldata3 = scv.read('/data3/yuhan/Project/Neuron/scMultiome/scMultiome_library/9_GW11B3_FU_scMultiome_ATAC_GEX_B1/velocyto/9_GW11B3_FU_scMultiome_ATAC_GEX_B1_velocyto/gex_possorted_bam_OLXOC.loom', cache=True)

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = ['GW20_' + bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = ['GW15_' + bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = ['GW11_' + bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata3.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()

# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color=['Combine_CellType'], legend_loc='on data', legend_fontsize=10, legend_fontweight='normal', size=20, frameon=False, title='', save='UMAP.Combine_CellType.check2.pdf')

# save dataset as anndata format
# adata.write('AllCells_addingLoom.h5ad')
# reload dataset
# adata = sc.read_h5ad('AllCells_addingLoom.h5ad')


###############################################################################################################################################################
# run.scVelo_step4.py
#!/usr/bin/env python
# coding=utf-8

###############################################################################################################################################################
conda activate scVelo_update
ipython

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False, figsize=(7,7))
cr.settings.verbosity = 2

###############################################################################################################################################################
# reload dataset
# adata = sc.read_h5ad('AllCells_addingLoom.h5ad')

# plot umap to check
sc.pl.umap(adata, color=['Combine_CellType'], legend_loc='on data', legend_fontsize=10, legend_fontweight='normal', size=20, frameon=False, title='', save='UMAP.Combine_CellType.check3.pdf')

###############################################################################################################################################################
# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
# Visualize velocity fields
# scv.pl.velocity_embedding(adata, basis='X_umap', frameon=False, save='UMAP.Combine_CellType.embedding.pdf')
adata.uns['Combine_CellType_colors'] = np.array(["#00441B","#C7E9C0","#A1D99B","#4575B4","#ABD9E9","#FDAE61","#F46D43","#D73027","#A50026","#FEE090","#41AB5D","#238B45","#006D2C","#313695","#FFFFBF","#FFFFBF","#FFFFBF","#74ADD1"])
scv.pl.velocity_embedding_grid(adata, basis='X_umap', color='Combine_CellType', save='UMAP.Combine_CellType.embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='X_umap', color='Combine_CellType', size=30, fontsize=0.5, save='UMAP.Combine_CellType.embedding_stream.pdf', title='')

# plot velocity of a selected gene
# scv.pl.velocity(adata, var_names=['EOMES'], color='Combine_CellType',save='EOMES.pdf')

# identify highly dynamic genes
# compute a measure of coherence among neighboring cells in terms of velocity
# perform pseudotime inference
# identify predicted ancestors of individual cells Using the pseudotime trajectory
# orient the directionality of partition-based graph abstractions (PAGA) Using the pseudotime trajectory
scv.tl.rank_velocity_genes(adata, groupby='Combine_CellType', min_corr=.3)
# df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
# df.head()
# kwargs = dict(frameon=False, size=10, linewidth=1.5,
#               add_outline='nIPC, RG')
# scv.pl.scatter(adata, df['nIPC'][:3], ylabel='nIPC', frameon=False, color='Combine_CellType', size=10, linewidth=1.5,save='nIPC.pdf')
# scv.pl.scatter(adata, df['RG'][:3], ylabel='RG', frameon=False, color='Combine_CellType', size=10, linewidth=1.5,save='RG.pdf')

scv.tl.velocity_confidence(adata)
# keys = 'velocity_length', 'velocity_confidence'
# scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save='p1.pdf')
# scv.pl.velocity_graph(adata, threshold=.1, color='Combine_CellType',save='p2.pdf')
# x, y = scv.utils.get_cell_transitions(adata, basis='X_umap', starting_cell=70)
# scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False,save='p3.pdf')
# scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax,save='p4.pdf')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', size=20, cmap='Purples', save='UMAP.Combine_CellType.pseudotime.pdf', title='')

###############################################################################################################################################################
# Analyzing a specific cell population
# "RG1","RG2","RG3","RG4","RG5","nIPC_IN","IN","CGE_IN1","MGE_IN","CGE_IN2","RG_Cyc","NPC",
# "nIPC_Cyc","nIPC_GluN","GluN1","GluN2","GluN3","GluN4","GluN5","GluN6","GluN1_SATB2","GluN2_SATB2","GluN3_SATB2","GluN4_SATB2",
# "SP","Oligo","MG","Peric"
# ,"GluN1_SATB2","GluN2_SATB2","GluN3_SATB2","GluN4_SATB2"

cur_celltypes = ["nIPC_Cyc","nIPC_GluN","GluN1","GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8","GluN1_SATB2","GluN2_SATB2","GluN3_SATB2","GluN4_SATB2"]
# save dataset as anndata format
# adata_subset.write('GluN.h5ad')
# reload dataset
# adata_subset = sc.read_h5ad('GluN.h5ad')

cur_celltypes = ["nIPC_GluN","GluN1","GluN2","GluN3","GluN4","GluN5","GluN6","GluN7","GluN8","GluN1_SATB2","GluN2_SATB2","GluN3_SATB2","GluN4_SATB2"]
# save dataset as anndata format
# adata_subset.write('GluN_subtype1.h5ad')
# reload dataset
# adata_subset = sc.read_h5ad('GluN_subtype1.h5ad')

# cur_celltypes = ["nIPC_Cyc","nIPC_GluN","GluN1_SATB2","GluN2_SATB2","GluN3_SATB2","GluN4_SATB2"]
# save dataset as anndata format
# adata_subset.write('GluN_subtype2.h5ad')
# reload dataset
# adata_subset = sc.read_h5ad('GluN_subtype2.h5ad')

adata_subset = adata[adata.obs['Combine_CellType'].isin(cur_celltypes)]
sc.pl.umap(adata_subset, color=['Combine_CellType'], legend_loc='on data', legend_fontsize=10, legend_fontweight='normal', size=20, frameon=False, title='', save='UMAP.Combine_CellType.check4.pdf')

# pre-process
scv.pp.filter_and_normalize(adata_subset,min_shared_counts=10,n_top_genes=1000)
scv.pp.moments(adata_subset,n_pcs=30,n_neighbors=30)

# compute velocity
scv.tl.velocity(adata_subset, mode='stochastic')
scv.tl.velocity_graph(adata_subset)
scv.pl.velocity_embedding_grid(adata_subset, basis='X_umap', color='Combine_CellType', save='GluN.UMAP.Combine_CellType.embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata_subset, basis='X_umap', color='Combine_CellType', size=30, fontsize=0.5, save='GluN.UMAP.Combine_CellType.embedding_stream.pdf', title='')

## Identify important genes
# scv.tl.rank_velocity_genes(adata_subset, groupby='Combine_CellType', min_corr=.3)
# scv.tl.velocity_confidence(adata_subset)

## Velocities in cycling progenitors
# scv.tl.score_genes_cell_cycle(adata_subset)
# scv.pl.scatter(adata_subset, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

# # Computes terminal states
# scv.tl.terminal_states(adata_subset)
# scv.pl.scatter(adata_subset, color=['root_cells', 'end_points'])

# pseudotime
scv.tl.velocity_pseudotime(adata_subset)
scv.pl.scatter(adata_subset, color='velocity_pseudotime', size=30, cmap='Oranges', save='GluN.UMAP.Combine_CellType.pseudotime.pdf', title='')

# PAGA
adata_subset.uns['neighbors']['distances'] = adata_subset.obsp['distances']
adata_subset.uns['neighbors']['connectivities'] = adata_subset.obsp['connectivities']
scv.tl.paga(adata_subset, groups='Combine_CellType')
scv.pl.paga(adata_subset, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5,figsize=(7,7), save='GluN.UMAP.Combine_CellType.PAGA.pdf', title='')
