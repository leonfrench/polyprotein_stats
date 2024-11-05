import gc
import time
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import functools
import numpy as np
import os, random
from datetime import datetime
from sklearn.metrics import f1_score, average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from scipy import stats as scistats
from statsmodels.stats import multitest
import bokeh.io
import bokeh.plotting

embedding_file_path_processed = os.path.join(os.path.dirname(__file__), 'data', 'processed')

#cache the file loading to speed things up
@st.cache_data
def get_file_with_cache(filename):
    df = pd.read_csv(os.path.join(embedding_file_path_processed, filename))
    return df

#this is split in two to allow easy deployment from github
@st.cache_data
def get_split_embeddings():
    dfA = get_file_with_cache("geneformer_gene_embeddings.csv.zip")
    #dfB = get_file_with_cache("gene_symbol_summarized_prottrans_t5_xl_u50.2.csv.zip")
    full = dfA
    return full 

#get download for predicting on everything
def get_download_button(X, y, all_embeddings, name):
    #create a genome-wide ranking by training on the target genes
    model = LogisticRegression()
    model.fit(X, y)
    # Extract predictions from fitted model and add gene symbol
    preds = model.predict(X)
    probas = pd.DataFrame(model.predict_proba(X), columns=model.classes_)
    probas['classification_target'] = all_embeddings['classification_target'].tolist()
    probas['gene_symbol'] = all_embeddings['gene_symbol'].tolist()
    probas = probas.sort_values(True, ascending=False)
    probas = probas.drop([False], axis=1)
    probas = probas.rename(columns = {True : "pvalue_for_target_set"})
    cols = probas.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    probas = probas[cols]
    st.download_button(
        "Download classification predictions for all proteins using " + name,
        probas.to_csv(index=False).encode('utf-8'),
        name + "_genome_wide_predictions.csv",
        "text/csv",
        key=name + '-download-csv'
    )

#copies are needed because it gets modified - helps with cacheing
all_embeddings = get_split_embeddings().copy()
embedding_UMAP = get_file_with_cache("gene_symbol_summarized_UMAP_Geneformer.csv").copy()
proportions = get_file_with_cache("gene_symbol_summarized_proportions.csv").copy()


st.sidebar.write("""### Geneformer gene embedding probe tool by Leon French
Embeddings for the human genes are from the Geneformer model weights from the [Hugging Face](https://huggingface.co/ctheodoris/Geneformer) repository. Full details on Geneformer by Theodoris et al. is available at [Nature](https://www.nature.com/articles/s41586-023-06139-9). The default gene list is from the [Lindbohm et al.](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12419) study of cognitive decline and risk of dementia.


""")


######################
# Input Text Boxes
######################

sequence_input = "CDCP1\nCHRDL1\nGDF15\nGM2A\nIGFBP7\nMARCKSL1\nNPPB\nPLA2G2A\nRNASE6\nSIGLEC7\nSVEP1\nTIMP4\nTNFRSF11B\nTREM2\nWFDC2"

target_genes = st.sidebar.text_area("Target genes input (gene symbol per line)", sequence_input, height=100)

#background_genes = st.sidebar.text_area("Background genes - leave blank for all protein coding genes", "SMAD1\nFKBP5\nMT-CO3\nARHGEF3\nNEAT1\nMT-ND4\nMT-ND3\nTENM4\nHSPA1A\nHMGB1\nCADPS\nCLIC4\nMT-ATP6\nOSBPL1A\nNRXN3\nTPST1\nCD44\nATP9A\nHSPB1\nSLC7A11\nPTGES3", height=100)
background_genes = st.sidebar.text_area("Background genes - leave blank for all protein coding genes", "", height=100)


st.sidebar.write("""
#### Motivation

Given a set of proteins, this tool seeks to answer these questions:

* is there a enrichment of amino acids composition? 
* can amino acid composition, plus sequence length discriminate the input proteins from the rest of the proteome?
* can the Geneformer weights, without finetuning discriminate the input genes from the rest of the proteome?

Source code is on [github](https://github.com/leonfrench/polyprotein_stats)
""")

target_genes = target_genes.splitlines()

background_matrices_genes = set(all_embeddings['gene_symbol']).intersection(proportions['gene_symbol'])

if (background_genes == ""):
  background_genes = background_matrices_genes
else:
  background_genes = background_genes.splitlines()


target_genes = set(target_genes)
target_genes_found =  target_genes.intersection(background_matrices_genes)

background_genes = set(background_genes)
background_genes_found =  background_genes.intersection(background_matrices_genes)


#could be a single dataframe
all_embeddings['classification_target'] = all_embeddings['gene_symbol'].isin(target_genes_found)
embedding_UMAP['classification_target'] = embedding_UMAP['gene_symbol'].isin(target_genes_found)
proportions['classification_target'] = proportions['gene_symbol'].isin(target_genes_found)
  
all_embeddings = all_embeddings[all_embeddings['gene_symbol'].isin(background_genes_found)]
embedding_UMAP = embedding_UMAP[embedding_UMAP['gene_symbol'].isin(background_genes_found)]
proportions = proportions[proportions['gene_symbol'].isin(background_genes_found)]
#print(all_embeddings.shape)
#print(embedding_UMAP.shape)
#print(proportions.shape)


#ensure same order so the folds and targets line up
all_embeddings = all_embeddings.sort_values('gene_symbol')
proportions = proportions.sort_values('gene_symbol')


st.write("""
#### Amino acid composition enrichment tests

The figure and table below show results from testing if the given proteins are enriched for the proportions of amino acid pairs relative to the background proteins. 
This is not testing bigrams but instead the composition or proportions of specific residues. In the figure, color represents the area under the receiver operating curve, 
using the combined proportions of the amino acid pairs. The below table provides p-values with and without Bonferroni multiple test correction. Length is additionally added to the table to check if the input proteins are longer than the background set. 

""")


#residue proportion testing
residues = set(proportions.columns.values).difference(('gene_symbol','length', 'classification_target'))


aa_AUC = []
neg = proportions[proportions['classification_target'] == 0]
pos = proportions[proportions['classification_target'] == 1]
for residueA in residues:
  for residueB in residues:
    if residueA <= residueB: #avoid computing full matrix - just the triangle
      auc = roc_auc_score(proportions['classification_target'], proportions[residueA] + proportions[residueB])
      pvalue = scistats.mannwhitneyu((pos[residueA] + pos[residueB]).tolist(), (neg[residueA] + neg[residueB]).tolist()).pvalue

      aa_AUC.append({'residue A': residueA, 'residue B' : residueB, 'auc' : auc, 'pvalue' : pvalue})

aa_AUC_df = pd.DataFrame(aa_AUC)


aa_AUC_df_square = (aa_AUC_df.pivot(index=['residue B'],columns='residue A', values="auc")
         .sort_index(level=[1,0]))
    

sns.set_theme()

fig = plt.figure(figsize=(10, 7))
aa_AUC_df_square = aa_AUC_df_square.fillna(0.5)
ax = sns.heatmap(aa_AUC_df_square, center = 0.5, cmap = 'vlag')

st.pyplot(fig, clear_figure = True)
plt.close("all")

#tag on length AUC value, could just be printed
auc_for_length = roc_auc_score(proportions['classification_target'], proportions['length'])
neg = proportions[proportions['classification_target'] == 0]
pos = proportions[proportions['classification_target'] == 1]

p_for_length = scistats.mannwhitneyu(neg['length'].tolist(), pos['length'].tolist()).pvalue
aa_AUC.append({'residue A': 'length', 'residue B' : '', 'auc' : auc_for_length, 'pvalue' : p_for_length, 'pvalue_bonf': float("NaN")})
aa_AUC_df = pd.DataFrame(aa_AUC)
aa_AUC_df.reset_index(inplace=True, drop=True)
aa_AUC_df.index = [""] * len(aa_AUC_df)


aa_AUC_df['pvalue_bonf'] = multitest.multipletests(aa_AUC_df['pvalue'].tolist(), method="bonferroni")[1]

aa_AUC_df = aa_AUC_df.sort_values('auc', ascending=False)

st.write(aa_AUC_df.style.format({'auc' : "{:.2f}", "pvalue": "{:.2g}", "pvalue_bonf": "{:.2g}"}))

###classification 
#should be equal to proportions target - needs checking
best_predicted_genes = []

n_splits = 4

if len(target_genes) >= n_splits*2:
    y = all_embeddings['classification_target']
    
    X = all_embeddings.drop(['classification_target', 'gene_symbol'], axis = 1)
    X_proportions = proportions.drop(['classification_target', 'gene_symbol'], axis = 1)
    

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=1)
    
    auc_scores = []
    auprc_scores = []
    auc_scores_proportions = []
    
    with st.spinner('Please wait...'):
      for i, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        
          X_train = X.iloc[train_idx, :]
          y_train = y.iloc[train_idx]
          X_test = X.iloc[test_idx, :]
          y_test = y.iloc[test_idx]
          
          #st.write("fold:" + str(i))
          
          X_proportions_train = X_proportions.iloc[train_idx, :]
          X_proportions_test = X_proportions.iloc[test_idx, :]
          y_train = y.iloc[train_idx]
          y_test = y.iloc[test_idx]
      
          #st.write("fold after mem:" + str(i))
          
          model = LogisticRegression()
          #st.write("fold after init:" + str(i))
  
          gc.collect() 
  
          model.fit(X_train, y_train)
  
          #st.write("fold after fit:" + str(i))
  
          # Extract predictions from fitted model
          # probs for classes ordered in same manner as model.classes_
          probas = pd.DataFrame(model.predict_proba(X_test), columns=model.classes_)
  
          # Get metrics for each model
          auc = roc_auc_score(y_test, probas[True])
          auc_scores.append(auc)
          #track the gene with the max prediction for the true class
          best_predicted_genes.append(all_embeddings.iloc[test_idx[probas.idxmax()[True]],:]['gene_symbol'])
  
          #run the model again with proportions instead of embeddings
          model.fit(X_proportions_train, y_train)
          probas = pd.DataFrame(model.predict_proba(X_proportions_test), columns=model.classes_)
          auc = roc_auc_score(y_test, probas[True])
          auc_scores_proportions.append(auc)

          #area under precision recall curve
          auprc = average_precision_score(y_test, probas[True])
          auprc_scores.append(auprc)
        
        
    best_predicted_genes = set(best_predicted_genes)
    top_predicted_hits = best_predicted_genes.intersection(target_genes)
    best_predicted_genes = best_predicted_genes.difference(target_genes)
    
    measures = {
                'number_of_input_genes': len(target_genes),
                'number_of_used_genes': len(target_genes_found),
                'number_of_input_background_genes': len(background_genes),
                'number_of_used_background_genes': len(background_genes_found),
                'area_under_precision_recall_curve_embeddings': np.mean(auprc_scores), 
                'AUC': np.mean(auc_scores),
                'AUC standard dev': np.std(auc_scores),
                'gain over proportions' : np.mean(auc_scores)-np.mean(auc_scores_proportions),
                'AUC_p_value versus 0.5' : scistats.ttest_1samp(auc_scores, 0.5).pvalue, #two sided p-value testing the AUC values against 0.5 expecation
                'AUC proportions': np.mean(auc_scores_proportions),
                'AUC proportions p_value versus 0.5' : scistats.ttest_1samp(auc_scores_proportions, 0.5).pvalue, #two sided p-value testing the AUC values against 0.5 expecation
                'AUC proportions vrs embeddings pvalue' : scistats.ttest_rel(auc_scores, auc_scores_proportions).pvalue,
                'top_false_positives_in_folds' : best_predicted_genes,
                'top_true_positives_in_folds' : top_predicted_hits}

    st.write("""#### Classification results

To test if the input proteins can be discriminated from the background proteins based on their residue proportions (plus length) and 
learned embeddings. This data derived from only sequences was used to train and test a logistic regression model (4 fold cross-validation, 
L2 loss, sklearn default parameters) that attempts to classify proteins as belonging to the input set. Given that the input genes are probably fewer than the background genes, we again report the AUC statistic.

""")
    #st.write('There are  ' + str(X['A']) + ' adenine (A)')

    st.markdown(f"Using proportions and length alone, the average AUC is **{measures['AUC proportions']:.2f}**.")

    st.markdown(f"Using the learned embeddings (see sidepanel for details), the average AUC is **{measures['AUC']:.2f}**.")
    
    st.markdown(f"Testing if the embedding based AUC values for the {n_splits} folds deviate from the expected 0.5 reveals a p-value of **{measures['AUC_p_value versus 0.5']:.2g}**.")    

    st.markdown(f"The difference between the AUC values from the complex embeddings and simple proportions is {measures['gain over proportions']:.2g} (p-value = **{measures['AUC proportions vrs embeddings pvalue']:.2g}**, paired t-test).")    
    
    st.write("More statistics from the classification tests are in the below dictionary:")
    st.write(measures)
    
    get_download_button(X, y, all_embeddings, "embeddings")
    #input, X, y, n_jobs, all_embeddings, name
    get_download_button(X_proportions, y, all_embeddings, "proportions")

    
else:
	st.write("#### Classification results")
	st.write("Too few genes to run classification task - skipping")




st.write("""#### Embedding visualization

The below plot shows the UMAP visualization of the embeddings of the input (red, target=true on hover) and background proteins (grey). The blue dots mark top predicted 
proteins for each fold of the classification model (also top_predicted_hits_in_folds in the above dictionary). These predicted genes may help understand what the classifier is learning from the input proteins. 
You can zoom and pan using the sidebar buttons. 

""")

## Bokeh Plot for UMAP

x='UMAP 0'
y='UMAP 1'
p = bokeh.plotting.figure(
    width=720,
    height=480,
    x_axis_label=x,
    y_axis_label=y,
    active_scroll="wheel_zoom",
    tooltips=[
        ("Gene", "@{gene_symbol}"),('Target', '@{classification_target}')
    ],
)
p.circle(
    source=embedding_UMAP[embedding_UMAP['classification_target'] == False], 
    x=x, y=y,
    fill_color='lightgrey', line_color = 'lightgrey'
)
#top 5 genes from the predictors, one for each fold
p.circle(
    source=embedding_UMAP[embedding_UMAP['gene_symbol'].isin(best_predicted_genes)], 
    x=x, y=y,
    fill_color='blue', line_color = 'blue'
)
p.circle(
    source=embedding_UMAP[embedding_UMAP['classification_target'] == True], 
    x=x, y=y,
    fill_color='red', line_color = 'red'
)
st.bokeh_chart(p, use_container_width=True)



