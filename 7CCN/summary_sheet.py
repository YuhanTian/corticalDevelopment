import pandas as pd
import numpy as np

import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#create output figure directory
model_name = 'test'
output_file_path = "../outputs/" + model_name + "/motifs/"
directory = os.path.dirname(output_file_path)
if not os.path.exists(directory):
    print("Creating directory %s" % output_file_path)
    os.makedirs(directory)
else:
     print("Directory %s exists" % output_file_path)

#files
meme_file = output_file_path + "filter_motifs_pwm.meme"
infl_file = output_file_path + "ave_filt_infl.npy"
cell_infl_file = output_file_path + "infl.npy"
nseqs_file = output_file_path + "nseqs.npy"

tomtom_file = output_file_path + "tomtom/tomtom.tsv"
  
#load cell type names 
col_names = ["RG_Cyc","RG1","RG2","RG3","RG4",
             "MigN","superficialCPN","Layer4CPN",
             "deepCPN","CThPN","nIPC1","nIPC2","nIPC3",
             "ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN"]
#filter names
rows = ['filter' + str(i) for i in range(300)]

#tomtom results
#load tomtom results
tomtom = pd.read_csv(tomtom_file, sep='\t')
tomtom = tomtom[:-3]

#read in motif TF information
#files downloaded from the CisBP database 
TF_info = pd.read_csv('../data/TF_Information_all_motifs.txt', sep='\t')

#leave only rows with q-value <0.05
tomtom = tomtom.loc[tomtom['q-value']<0.05]

#merge results with TF information 
tomtom = tomtom.merge(TF_info[['Motif_ID', 'TF_Name', 'Family_Name', 'TF_Status', 'Motif_Type']], how='left', left_on='Target_ID', right_on='Motif_ID')
tomtom = tomtom.drop(columns=['Motif_ID']) 

# get top 5 transcription factor matches for each motif
temp = tomtom[['Target_ID', 'TF_Name', 'Family_Name', 'TF_Status', 'Motif_Type']].drop_duplicates()
temp = temp.groupby(['Target_ID']).agg(lambda x : ', '.join(x.astype(str)))

# combine matrices
tomtom = tomtom.drop(['TF_Name', 'Family_Name', 'TF_Status', 'Motif_Type'], axis=1)
tomtom = tomtom.merge(temp, how='left', left_on=['Target_ID'], right_on=['Target_ID']).drop_duplicates()

tomtom['log_qval'] = -np.log(tomtom['q-value'])
tomtom['log_qval'].fillna(0, inplace=True)

#load influence score
infl = np.load(infl_file)

rows = ['filter' + str(i) for i in range(infl.shape[0])]
infl = pd.DataFrame(infl, columns=["Influence"])
infl.index = rows 
tomtom = tomtom.merge(infl, how='outer', left_on='Query_ID',right_index=True)

#load cell-wise influence for individual filters
cell_infl = np.load(cell_infl_file)
cell_infl = pd.DataFrame(cell_infl, index=rows, columns = col_names)
tomtom = tomtom.merge(cell_infl, how='outer', left_on='Query_ID', right_index=True)

# number of sequences per filter:
nseqs = np.load(nseqs_file)
nseqs = pd.DataFrame(nseqs, columns=['nseqs'])
nseqs.index = rows

#add to data frame
tomtom = tomtom.merge(nseqs, how='outer', left_on='Query_ID', right_index=True)


#read in meme file
with open(meme_file) as fp:
    line = fp.readline()
    motifs=[]
    motif_names=[]
    while line:
        #determine length of next motif
        if line.split(" ")[0]=='MOTIF':
            #add motif number to separate array
            motif_names.append(line.split(" ")[1])
            
            #get length of motif
            line2=fp.readline().split(" ")
            motif_length = int(float(line2[5]))
            
            #read in motif 
            current_motif=np.zeros((19, 4))
            for i in range(motif_length):
                current_motif[i,:] = fp.readline().split("\t")
            
            motifs.append(current_motif)

        line = fp.readline()
        
    motifs = np.stack(motifs)
    motif_names = np.stack(motif_names)

fp.close()

#set background frequencies of nucleotides
bckgrnd = [0.25, 0.25, 0.25, 0.25]

#compute information content of each motif
info_content = []
position_ic = []
for i in range(motifs.shape[0]):
    length = motifs[i,:,:].shape[0]
    position_wise_ic = np.subtract(np.sum(np.multiply(motifs[i,:,:],np.log2(motifs[i,:,:] + 0.00000000001)), axis=1),np.sum(np.multiply(bckgrnd,np.log2(bckgrnd))))
                                        
    position_ic.append(position_wise_ic)
    ic = np.sum(position_wise_ic, axis=0)
    info_content.append(ic)
    
info_content = np.stack(info_content)
position_ic = np.stack(position_ic)

#length of motif with high info content
n_info = np.sum(position_ic>0.2, axis=1)

#"length of motif", i.e. difference between first and last informative base
ic_idx = pd.DataFrame(np.argwhere(position_ic>0.2), columns=['row', 'idx']).groupby('row')['idx'].apply(list)
motif_length = []
for row in ic_idx:
    motif_length.append(np.max(row)-np.min(row) + 1)

motif_length = np.stack(motif_length)

#create pandas data frame:
info_content_df = pd.DataFrame(data=[motif_names, info_content, n_info, pd.to_numeric(motif_length)]).T
info_content_df.columns=['Filter', 'IC', 'Num_Informative_Bases', 'Motif_Length']

tomtom = tomtom.merge(info_content_df, how='left', left_on='Query_ID', right_on='Filter')
tomtom = tomtom.drop(columns=['Filter'])
tomtom['IC'] = pd.to_numeric(tomtom['IC'])
tomtom['Num_Informative_Bases'] = pd.to_numeric(tomtom['Num_Informative_Bases'])

#compute consensus sequence
for filter in tomtom[pd.isnull(tomtom['Query_consensus'])]['Query_ID']:
    sequence = ''
    indx = np.argwhere(motif_names == filter)
    pwm = motifs[indx,:,:].squeeze()
    for j in range(pwm.shape[0]):
        letter = np.argmax(pwm[j, :])
        if(letter==0):
            sequence = sequence + 'A'
        if(letter==1):
            sequence = sequence + 'C'
        if(letter==2):
            sequence = sequence + 'G'
        if(letter==3):
            sequence = sequence + 'T'
            
    tomtom.loc[(tomtom['Query_ID']==filter), 'Query_consensus'] = sequence


tomtom.to_csv("../outputs/test/motifs/motif_summary.csv")

################################################################
tomtom = pd.read_csv("../outputs/test/motifs/motif_summary.csv", sep=',', header=0)
filter_all_information = tomtom[['Query_ID','Influence','nseqs','IC','Num_Informative_Bases','Motif_Length']]
filter_all_information = filter_all_information.drop_duplicates()

filter_enrich_counts = np.load('../outputs/test/motifs/filter_enrich_counts.npy',allow_pickle=True)
filter_enrich_counts = pd.DataFrame(filter_enrich_counts, columns=['Query_ID','Fold10_counts'])

filter_all_information = filter_all_information.merge(filter_enrich_counts, how='outer', left_on='Query_ID', right_on='Query_ID')
filter_all_information.fillna(0, inplace=True) 
filter_all_information['Fold10_counts'] = filter_all_information['Fold10_counts']+1
filter_all_information['Influence'] = np.log10(filter_all_information['Influence'])

import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['xtick.labelsize'] = 4.5
mpl.rcParams['ytick.labelsize'] = 4.5
mpl.rcParams['axes.labelsize'] = 6
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['legend.title_fontsize'] = 6

palette = ["#1A328A","#4E61AD","#489CB4","#88D0A8","#CAE99D","#F7FBAF","#FEE99A","#FCB76C","#F47646","#D5414F","#9D0142"]

## IC_vs_Infl_counts
plt.figure(figsize=(1.8, 1.8))
fig = sns.scatterplot(data=filter_all_information, x='IC', y='Influence', 
                      hue='Fold10_counts', size=6, palette=palette)
fig.set_xlim(0,17)
fig.set_ylim(-6,-1)
fig.xaxis.label.set_size(6)
fig.yaxis.label.set_size(6)
plt.xticks(fontsize=4.5)
plt.yticks(fontsize=4.5)
plt.legend(bbox_to_anchor=(1.3, 1), markerscale=0.4, scatterpoints=1)
plt.savefig(output_file_path + "IC_vs_Infl_counts.pdf")

## Infl_vs_nseqs_counts
tomtom_motif_v1 = pd.read_csv("../outputs/test/motifs/motif_summary_filter_v1.csv", sep=',', header=0)
filter_all_information_v1 = tomtom_motif_v1[['Query_ID','TF_Name','Influence','nseqs','IC','Num_Informative_Bases','Motif_Length']]
Query_ID = list(filter_all_information_v1['Query_ID'])
for idx, ele in enumerate(Query_ID):
    Query_ID[idx] = ele.replace('filter','')
filter_all_information_v1['Query_ID_anno'] = pd.Series(Query_ID)
filter_all_information_v1.fillna('', inplace=True)
filter_all_information_v1['Query_ID_anno'] = filter_all_information_v1['Query_ID_anno'] + ' ' + filter_all_information_v1['TF_Name']

filter_all_information_v1 = filter_all_information_v1.merge(filter_enrich_counts, how='outer', left_on='Query_ID', right_on='Query_ID')
filter_all_information_v1.fillna(0, inplace=True) 
filter_all_information_v1['Fold10_counts'] = filter_all_information_v1['Fold10_counts']+1
filter_all_information_v1['Influence'] = np.log10(filter_all_information_v1['Influence'])
filter_all_information_v1['nseqs'] = np.log10(filter_all_information_v1['nseqs'])

filter_all_information_v1_highcounts = filter_all_information_v1.loc[filter_all_information_v1['Fold10_counts'] >= 6,:]
filter_all_information_v1_highcounts = filter_all_information_v1_highcounts.reset_index(drop=True)

palette_2 = ["#F7FBAF","#FEE99A","#FCB76C","#F47646","#D5414F","#9D0142"]

plt.figure(figsize=(3.3, 3.3))
fig = sns.scatterplot(data=filter_all_information_v1_highcounts,
                      x='Influence', y='nseqs', size=6, 
                      hue = 'Fold10_counts', palette=palette_2)

fig.xaxis.label.set_size(6)
fig.yaxis.label.set_size(6)
plt.xticks(fontsize=4.5)
plt.yticks(fontsize=4.5)
plt.legend(bbox_to_anchor=(1.3, 1), markerscale=0.4, scatterpoints=1)

for i in range(filter_all_information_v1_highcounts.shape[0]):
    plt.text(x=filter_all_information_v1_highcounts.Influence[i]-0.1,
             y=filter_all_information_v1_highcounts.nseqs[i]-0.1,
             s=filter_all_information_v1_highcounts.Query_ID_anno[i],
             fontdict=dict(color='black',size=6))
 
plt.savefig(output_file_path + "Infl_vs_nseqs_counts.pdf")

# Infl in cell type
predictions = np.load(output_file_path + "predictions.npy")
pred_full_model2 = np.load(output_file_path + "pred_full_model2.npy")
n_combos = predictions.shape[1]
pred_full_model2 = np.expand_dims(pred_full_model2, 1)
pred_full_model2 = np.repeat(pred_full_model2, n_combos, axis=1)

diff_pred = predictions - pred_full_model2
infl_filter_cell_type = np.mean(diff_pred, axis=0).squeeze()
highcounts_index = filter_all_information_v1_highcounts['Query_ID']
for idx, ele in enumerate(highcounts_index):
    highcounts_index[idx] = ele.replace('filter','')
infl_repro_filter_cell_type = pd.DataFrame(infl_filter_cell_type[highcounts_index.astype(str).astype(int)])
infl_repro_filter_cell_type.index = filter_all_information_v1_highcounts['Query_ID_anno']
col_names = ["RG_Cyc","RG1","RG2","RG3","RG4",
             "MigN","superficialCPN","Layer4CPN",
             "deepCPN","CThPN","nIPC1","nIPC2","nIPC3",
             "ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN"]
infl_repro_filter_cell_type.columns = col_names

from matplotlib.colors import LinearSegmentedColormap
clist = ['#228544','white','#CB2D27']
newcmp = LinearSegmentedColormap.from_list('chaos',clist,N=256)

plt.figure(figsize=(9.3, 2))
sns.heatmap(infl_repro_filter_cell_type.T, cmap=newcmp)
plt.title("Influence of filters on each cell type prediction")
plt.xlabel("Filter")
plt.ylabel("Cell Type")
plt.savefig(output_file_path + "filt_infl_celltype_heatmap.pdf")

# plt.yticks(np.arange(len(col_names)), col_names, fontsize=3.0)
# plt.xticks(fontsize=3.0)

# top 500 OCRs for each filter
import heapq
activations = np.load(output_file_path + "activations.npy")
y2 = np.load(output_file_path + "y2.npy")
OCR_matrix = np.load(output_file_path + "OCR_matrix.npy")

filter_cell_type = np.zeros((300, y2.shape[1]))
for filter_index in np.arange(300):
    activations_OCR = list(np.sum(activations[:,filter_index,:], axis=1))
    activations_OCR_top_index = list(map(activations_OCR.index, heapq.nlargest(500,activations_OCR)))
    y2_filter = np.mean(y2[activations_OCR_top_index,:], axis=0)
    filter_cell_type[filter_index,:] = y2_filter
    
filter_TF = pd.read_csv(output_file_path + "filter_TF_final_v2.txt", sep=r'\s+', header=0)
filter_TF = filter_TF[((~pd.isna(filter_TF['TF'])) & (~(filter_TF['TF']=='BORCS8-MEF2B')))]
filter_TF = filter_TF.reset_index(drop=True)
TF_RNA = pd.read_csv(output_file_path + "TF_RNA_matrix_v2.txt", sep=r'\s+', header=0)
TF_RNA = TF_RNA[['TF',"RG_Cyc","RG1","RG2","RG3","RG4",
             "MigN","superficialCPN","Layer4CPN",
             "deepCPN","CThPN","nIPC1","nIPC2","nIPC3",
             "ImmatureN1","ImmatureN2","SCPN","SCPN_CSMN","SCPN_CTPN"]]

filter_TF['cor'] = 0 
for TF_index in np.arange(filter_TF.shape[0]):
    filter = filter_TF.loc[TF_index,'filter']
    TF = filter_TF.loc[TF_index,'TF']
    TF_exp = TF_RNA.loc[TF_RNA['TF'] == TF,:]
    TF_exp = np.array(TF_exp.iloc[0,1:19], dtype='float64')
    filter_exp = filter_cell_type[int(filter.replace('filter','')),:]
    cor = np.corrcoef(TF_exp, filter_exp)[0,1]
    filter_TF.loc[TF_index,'cor'] = cor 

filter_TF_plot = filter_TF[~pd.isna(filter_TF['cor'])]
filter_TF_plot = filter_TF_plot.drop_duplicates(subset=['filter','TF'],keep='first',inplace=False)
Fold10counts = filter_all_information_v1[['Query_ID','Fold10_counts']]
filter_TF_plot = filter_TF_plot.merge(Fold10counts, how='left', left_on='filter', right_on='Query_ID')

palette = ["#1A328A","#4E61AD","#489CB4","#88D0A8","#CAE99D","#F7FBAF","#FEE99A","#FCB76C","#F47646","#D5414F","#9D0142"]
plt.figure(figsize=(7, 3))
fig = sns.scatterplot(data=filter_TF_plot,
                      x='filter', y='cor', size=6,
                      hue='Fold10_counts', palette=palette)
fig.xaxis.label.set_size(8)
fig.yaxis.label.set_size(8)
plt.xticks(rotation='vertical', fontsize=6)
plt.yticks(fontsize=6) 
plt.legend(bbox_to_anchor=(1, 1), markerscale=0.4, scatterpoints=1)
plt.savefig(output_file_path + "filter_cor_gene_express.pdf")

filter_TF_plot_sort = filter_TF_plot.groupby('filter').apply(lambda x: x.sort_values(by=["cor"], ascending = False)).reset_index(drop=True)
filter_TF_plot_sort.to_csv(output_file_path + "filter_cor_gene_express.csv")







