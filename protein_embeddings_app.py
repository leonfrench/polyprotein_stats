import gc
import streamlit as st
import pandas as pd
import numpy as np
import os
import re
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from scipy import stats as scistats

st.set_page_config(layout="wide")
embedding_file_path_processed = os.path.join(os.path.dirname(__file__), 'data', 'processed')
MIN_POSITIVE_FLOAT = np.nextafter(0.0, 1.0)

#cache the file loading to speed things up
@st.cache_data
def get_file_with_cache(filename):
    df = pd.read_csv(os.path.join(embedding_file_path_processed, filename))
    return df

@st.cache_data
def get_paxdb_with_cache(filename):
    df = pd.read_csv(
        os.path.join(embedding_file_path_processed, filename),
        sep="\t",
        comment="#",
        header=None,
        names=['gene_symbol', 'string_external_id', 'abundance']
    )
    return df

def get_auc_and_pvalue(df, value_col):
    y = df['classification_target']
    if y.nunique() < 2:
        return float("nan"), float("nan"), int(y.sum()), int(len(y))
    auc = roc_auc_score(y, df[value_col])
    neg_values = df.loc[y == 0, value_col].dropna()
    pos_values = df.loc[y == 1, value_col].dropna()
    if len(neg_values) < 2 or len(pos_values) < 2:
        return auc, float("nan"), int(y.sum()), int(len(y))
    pvalue = scistats.mannwhitneyu(
        neg_values.tolist(),
        pos_values.tolist(),
        alternative='two-sided',
        method = 'asymptotic'
    ).pvalue
    # scipy can underflow to 0 for extremely small p-values; display the smallest
    # positive float instead so the UI never shows an impossible p-value of 0.
    pvalue = max(float(pvalue), MIN_POSITIVE_FLOAT)
    return auc, pvalue, int(y.sum()), int(len(y))


def prepare_single_feature_table(df, source_col, feature_name=None):
    feature_name = feature_name or source_col
    feature_df = df[['gene_symbol', source_col]].copy()
    feature_df[source_col] = pd.to_numeric(feature_df[source_col], errors='coerce')
    feature_df = feature_df.dropna(subset=['gene_symbol', source_col])
    feature_df = feature_df.drop_duplicates(subset=['gene_symbol'], keep='first')
    return feature_df.rename(columns={source_col: feature_name})


def get_gene_set_enrichment_result(target_genes, hit_genes, background_genes):
    background_set = set(background_genes)
    target_set = set(target_genes).intersection(background_set)
    reference_set = background_set.difference(target_set)
    hit_set = set(hit_genes).intersection(background_set)

    target_hits = len(target_set.intersection(hit_set))
    reference_hits = len(reference_set.intersection(hit_set))
    if not target_set or not reference_set:
        return (
            float("nan"),
            'undetermined',
            target_hits,
            len(target_set),
            reference_hits,
            len(reference_set)
        )

    contingency_table = [
        [target_hits, len(target_set) - target_hits],
        [reference_hits, len(reference_set) - reference_hits]
    ]
    pvalue = scistats.fisher_exact(
        contingency_table,
        alternative='two-sided'
    ).pvalue
    pvalue = max(float(pvalue), MIN_POSITIVE_FLOAT)
    direction_comparison = (
        target_hits * len(reference_set)
        - reference_hits * len(target_set)
    )
    if direction_comparison > 0:
        direction = 'enrichment'
    elif direction_comparison < 0:
        direction = 'depletion'
    else:
        direction = 'no directional difference'
    return (
        pvalue,
        direction,
        target_hits,
        len(target_set),
        reference_hits,
        len(reference_set)
    )

#copies are needed because it gets modified - helps with cacheing
proportions = get_file_with_cache("gene_symbol_summarized_proportions.csv").copy()
gene_locations = get_file_with_cache("gene_symbol_locations.csv").copy()
avg_bulk_cpm = get_file_with_cache("average_bulk_CPM_91_experiments.csv").copy()
avg_bulk_cpm = avg_bulk_cpm.rename(columns={'hgnc_symbol': 'gene_symbol'})
de_prior_rank = get_file_with_cache("pnas.1802973116.sd02.updated.txt.csv").copy()
de_prior_rank = de_prior_rank.rename(columns={'gene_symbol': 'gene_symbol'})
multifunctional_go = get_file_with_cache("multifunctional.GO_2025.csv").copy()
multifunctional_go = multifunctional_go.rename(columns={'Gene': 'gene_symbol'})
shared_dim1 = get_file_with_cache("GCCA_shared_embedding.center_and_normalize_False.dim_1.csv").copy()
gc_content = get_file_with_cache("GC_content.csv").copy()
inflammatome_rank = get_file_with_cache("Cort_et_al.inflammatome.all_constrasts.csv").copy()
inflammatome_rank['rank'] = pd.to_numeric(inflammatome_rank['rank'], errors='coerce')
inflammatome_rank = inflammatome_rank.dropna(subset=['gene_symbol', 'rank'])
inflammatome_rank = inflammatome_rank.sort_values('rank').drop_duplicates(subset=['gene_symbol'], keep='first')
# Lower rank means higher inflammatome signal; negate so AUROC direction is intuitive.
inflammatome_rank['inflammatome_score'] = -inflammatome_rank['rank']
homlof_populations = get_file_with_cache("Koch_et_al.homLoF.csv").copy()
if 'gene_symbol' not in shared_dim1.columns:
    shared_dim1 = shared_dim1.rename(columns={shared_dim1.columns[0]: 'gene_symbol'})
if 'dim1' not in shared_dim1.columns and len(shared_dim1.columns) > 1:
    shared_dim1 = shared_dim1.rename(columns={shared_dim1.columns[1]: 'dim1'})
if 'dim1' not in shared_dim1.columns and shared_dim1.shape[1] == 1:
    split_values = shared_dim1.iloc[:, 0].astype(str).str.split(r'[\s,]+', n=1, expand=True)
    if split_values.shape[1] > 1:
        shared_dim1 = pd.DataFrame(
            {
                'gene_symbol': split_values[0],
                'dim1': pd.to_numeric(split_values[1], errors='coerce')
            }
        )
paxdb_abundance = get_paxdb_with_cache("PaxDb_9606-WHOLE_ORGANISM-integrated.txt").copy()


st.sidebar.markdown("### Generic gene property tester")
st.sidebar.markdown(
    """
<style>
.sidebar-compact-list { font-size: 0.85rem; line-height: 1.2; margin-bottom: 0.6rem; }
.sidebar-compact-title { font-size: 1rem; margin-bottom: 0.2rem; }
.sidebar-compact-list ul { margin: 0.05rem 0 0.4rem 1.2rem; padding: 0; }
.sidebar-compact-list li { margin: 0; padding: 0; }
</style>
<div class="sidebar-compact-list">
  <div class="sidebar-compact-title">Given a set of input genes, this tool seeks to answer these questions:</div>
  <ul>
    <li>Does the gene set skew toward <strong>high average bulk expression</strong>?</li>
    <li>Does it favor genes with <strong>high expectation of differential expression</strong>?</li>
    <li>Are the genes more <strong>multifunctional</strong> than expected?</li>
    <li>Is there a detectable signal in a <strong>learned 1D representation</strong>?</li>
    <li>Do the corresponding proteins tend to have <strong>higher abundance</strong>?</li>
    <li>Are the genes <strong>more GC-rich</strong> than expected?</li>
    <li>Do they skew toward a <strong>high inflammatome signal</strong>?</li>
    <li>Are they enriched or depleted for <strong>observed human knockouts</strong>?</li>
    <li>Can <strong>amino acid composition + length</strong> separate the input proteins from the rest of the proteome?</li>
    <li>Can simple <strong>genomic location embeddings</strong> separate the input genes from the background?</li>
  </ul>
</div>
""",
    unsafe_allow_html=True
)


######################
# Input Text Boxes
######################

sequence_input = "CDCP1\nCHRDL1\nGDF15\nGM2A\nIGFBP7\nMARCKSL1\nNPPB\nPLA2G2A\nRNASE6\nSIGLEC7\nSVEP1\nTIMP4\nTNFRSF11B\nTREM2\nWFDC2"

target_genes = st.sidebar.text_area("Target gene list", sequence_input, height=100)

#background_genes = st.sidebar.text_area("Background genes - leave blank for all protein coding genes", "SMAD1\nFKBP5\nMT-CO3\nARHGEF3\nNEAT1\nMT-ND4\nMT-ND3\nTENM4\nHSPA1A\nHMGB1\nCADPS\nCLIC4\nMT-ATP6\nOSBPL1A\nNRXN3\nTPST1\nCD44\nATP9A\nHSPB1\nSLC7A11\nPTGES3", height=100)
background_genes = st.sidebar.text_area("Background genes (recommended)", "", height=100)

st.sidebar.markdown("Source code is on [github](https://github.com/leonfrench/polyprotein_stats/tree/generic_tester).")
st.sidebar.markdown(
    "The default gene list is from the [Lindbohm et al.](https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.12419) "
    "study of cognitive decline and risk of dementia."
)


target_genes = [gene for gene in re.split(r'[\s,;]+', target_genes.strip()) if gene]

background_matrices_genes = set(proportions['gene_symbol']).intersection(
    gene_locations['gene_symbol']
)


if (background_genes == ""):
  background_genes_input = None
else:
  background_genes_input = set(
      gene for gene in re.split(r'[\s,;]+', background_genes.strip()) if gene
  )


target_genes = set(target_genes)
target_genes_found =  target_genes.intersection(background_matrices_genes)
single_gene_input = len(target_genes) == 1

if background_genes_input is None:
  background_genes = background_matrices_genes
else:
  background_genes = background_genes_input

background_genes = set(background_genes)
background_genes_found =  background_genes.intersection(background_matrices_genes)


#could be a single dataframe
proportions['classification_target'] = proportions['gene_symbol'].isin(target_genes_found)
gene_locations['classification_target'] = gene_locations['gene_symbol'].isin(target_genes_found)

proportions_full = proportions.copy()
avg_bulk_cpm_full = avg_bulk_cpm.copy()
de_prior_rank_full = de_prior_rank.copy()
multifunctional_go_full = multifunctional_go.copy()
shared_dim1_full = shared_dim1.copy()
paxdb_abundance_full = paxdb_abundance.copy()
gc_content_full = gc_content.copy()
inflammatome_rank_full = inflammatome_rank.copy()
homlof_populations_full = homlof_populations.copy()
  
proportions = proportions[proportions['gene_symbol'].isin(background_genes_found)]
gene_locations = gene_locations[gene_locations['gene_symbol'].isin(background_genes_found)]

# print(proportions.shape)
# print(gene_locations.shape)


#ensure same order so the folds and targets line up
proportions = proportions.sort_values('gene_symbol')
gene_locations = gene_locations.sort_values('gene_symbol')
avg_bulk_cpm = avg_bulk_cpm.sort_values('gene_symbol')


st.write("""
#### Gene property enrichment tests

The table below shows whether the given gene set is biased toward or against mundane, basic, or non-specific signals. 

""")


if background_genes_input is None:
  table_proportions = proportions_full.copy()
  table_avg_bulk_cpm = avg_bulk_cpm_full.copy()
  table_de_prior_rank = de_prior_rank_full.copy()
  table_multifunctional_go = multifunctional_go_full.copy()
  table_shared_dim1 = shared_dim1_full.copy()
  table_paxdb_abundance = paxdb_abundance_full.copy()
  table_gc_content = gc_content_full.copy()
  table_inflammatome_rank = inflammatome_rank_full.copy()
else:
  table_proportions = proportions_full[proportions_full['gene_symbol'].isin(background_genes_input)].copy()
  table_avg_bulk_cpm = avg_bulk_cpm_full[avg_bulk_cpm_full['gene_symbol'].isin(background_genes_input)].copy()
  table_de_prior_rank = de_prior_rank_full[de_prior_rank_full['gene_symbol'].isin(background_genes_input)].copy()
  table_multifunctional_go = multifunctional_go_full[multifunctional_go_full['gene_symbol'].isin(background_genes_input)].copy()
  table_shared_dim1 = shared_dim1_full[shared_dim1_full['gene_symbol'].isin(background_genes_input)].copy()
  table_paxdb_abundance = paxdb_abundance_full[paxdb_abundance_full['gene_symbol'].isin(background_genes_input)].copy()
  table_gc_content = gc_content_full[gc_content_full['gene_symbol'].isin(background_genes_input)].copy()
  table_inflammatome_rank = inflammatome_rank_full[inflammatome_rank_full['gene_symbol'].isin(background_genes_input)].copy()

table_avg_bulk_cpm_protein = avg_bulk_cpm_full[
    avg_bulk_cpm_full['gene_symbol'].isin(proportions_full['gene_symbol'])
].copy()
if background_genes_input is not None:
  table_avg_bulk_cpm_protein = table_avg_bulk_cpm_protein[
      table_avg_bulk_cpm_protein['gene_symbol'].isin(background_genes_input)
  ].copy()

table_proportions['classification_target'] = table_proportions['gene_symbol'].isin(target_genes)
table_avg_bulk_cpm['classification_target'] = table_avg_bulk_cpm['gene_symbol'].isin(target_genes)
table_avg_bulk_cpm_protein['classification_target'] = table_avg_bulk_cpm_protein['gene_symbol'].isin(target_genes)
table_de_prior_rank['classification_target'] = table_de_prior_rank['gene_symbol'].isin(target_genes)
table_multifunctional_go['classification_target'] = table_multifunctional_go['gene_symbol'].isin(target_genes)
table_shared_dim1['classification_target'] = table_shared_dim1['gene_symbol'].isin(target_genes)
table_paxdb_abundance['classification_target'] = table_paxdb_abundance['gene_symbol'].isin(target_genes)
table_gc_content['classification_target'] = table_gc_content['gene_symbol'].isin(target_genes)
table_inflammatome_rank['classification_target'] = table_inflammatome_rank['gene_symbol'].isin(target_genes)

# Shared-dim1 is intentionally excluded for this combined ranking model.
ranking_feature_specs = [
    ('length', table_proportions, 'length'),
    ('avg_bulk_CPM', table_avg_bulk_cpm, 'avg_bulk_CPM'),
    ('DE_Prior_Rank', table_de_prior_rank, 'DE_Prior_Rank'),
    ('MF.score', table_multifunctional_go, 'MF.score'),
    ('abundance', table_paxdb_abundance, 'abundance'),
    ('gc_content', table_gc_content, 'gc_content'),
    ('inflammatome_score', table_inflammatome_rank, 'inflammatome_score'),
]

table_rankings_intersection = None
for feature_name, source_df, source_col in ranking_feature_specs:
    feature_table = prepare_single_feature_table(
        source_df,
        source_col=source_col,
        feature_name=feature_name
    )
    if table_rankings_intersection is None:
        table_rankings_intersection = feature_table
    else:
        table_rankings_intersection = table_rankings_intersection.merge(
            feature_table,
            on='gene_symbol',
            how='inner'
        )

table_rankings_intersection['classification_target'] = table_rankings_intersection['gene_symbol'].isin(target_genes)

#tag on length AUC value, could just be printed
auc_for_length, p_for_length, pos_length, total_length = get_auc_and_pvalue(table_proportions, 'length')
auc_for_cpm, p_for_cpm, pos_cpm, total_cpm = get_auc_and_pvalue(table_avg_bulk_cpm, 'avg_bulk_CPM')
auc_for_cpm_protein, p_for_cpm_protein, pos_cpm_protein, total_cpm_protein = get_auc_and_pvalue(
    table_avg_bulk_cpm_protein,
    'avg_bulk_CPM'
)
auc_for_de_prior, p_for_de_prior, pos_de_prior, total_de_prior = get_auc_and_pvalue(
    table_de_prior_rank,
    'DE_Prior_Rank'
)
auc_for_multifunctional, p_for_multifunctional, pos_multifunctional, total_multifunctional = get_auc_and_pvalue(
    table_multifunctional_go,
    'MF.score'
)
auc_for_shared_dim1, p_for_shared_dim1, pos_shared_dim1, total_shared_dim1 = get_auc_and_pvalue(
    table_shared_dim1,
    'dim1'
)
auc_for_paxdb, p_for_paxdb, pos_paxdb, total_paxdb = get_auc_and_pvalue(
    table_paxdb_abundance,
    'abundance'
)
auc_for_gc_content, p_for_gc_content, pos_gc_content, total_gc_content = get_auc_and_pvalue(
    table_gc_content,
    'gc_content'
)
auc_for_inflammatome, p_for_inflammatome, pos_inflammatome, total_inflammatome = get_auc_and_pvalue(
    table_inflammatome_rank,
    'inflammatome_score'
)
if background_genes_input is None:
    homlof_enrichment_background = set(proportions_full['gene_symbol'].dropna())
    homlof_background_label = 'protein-coding background'
else:
    homlof_enrichment_background = set(background_genes_input)
    homlof_background_label = 'user-provided background'

(
    homlof_fisher_pvalue,
    homlof_enrichment_direction,
    homlof_target_hits,
    homlof_target_total,
    homlof_reference_hits,
    homlof_reference_total
) = get_gene_set_enrichment_result(
    target_genes=target_genes,
    hit_genes=homlof_populations_full['gene_symbol'],
    background_genes=homlof_enrichment_background
)
aa_summary_df = pd.DataFrame(
    [
        {
            'Name': 'Protein coding length',
            'Source': '',
            'AUROC': auc_for_length,
            'pvalue': p_for_length,
            'Target genes used': pos_length,
            'Number total genes': total_length
        },
        {
            'Name': 'Average bulk gene expression',
            'Source': '<a href="https://cocoblast.ccbr.utoronto.ca/" target="_blank">CoCoBLAST</a>',
            'AUROC': auc_for_cpm,
            'pvalue': p_for_cpm,
            'Target genes used': pos_cpm,
            'Number total genes': total_cpm
        },
        {
            'Name': 'Average bulk gene expression within proteins',
            'Source': '<a href="https://cocoblast.ccbr.utoronto.ca/" target="_blank">CoCoBLAST</a>',
            'AUROC': auc_for_cpm_protein,
            'pvalue': p_for_cpm_protein,
            'Target genes used': pos_cpm_protein,
            'Number total genes': total_cpm_protein
        },
        {
            'Name': 'Differential expression prior',
            'Source': '<a href="https://pubmed.ncbi.nlm.nih.gov/30846554/" target="_blank">Crow et al.</a>',
            'AUROC': auc_for_de_prior,
            'pvalue': p_for_de_prior,
            'Target genes used': pos_de_prior,
            'Number total genes': total_de_prior
        },
        {
            'Name': 'GO multifunctional score',
            'Source': '<a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0017258" target="_blank">Gillis and Pavlidis</a>',
            'AUROC': auc_for_multifunctional,
            'pvalue': p_for_multifunctional,
            'Target genes used': pos_multifunctional,
            'Number total genes': total_multifunctional
        },
        {
            'Name': 'Shared 1D representation from large models',
            'Source': (
                '<a href="https://pubmed.ncbi.nlm.nih.gov/34232869/" target="_blank">ProtT5</a>'
                ', <a href="https://www.nature.com/articles/s41592-024-02201-0" target="_blank">scGPT</a>'
                ', <a href="https://www.biorxiv.org/content/10.1101/2024.10.10.617658v3" target="_blank">Orthrus</a>'
            ),
            'AUROC': auc_for_shared_dim1,
            'pvalue': p_for_shared_dim1,
            'Target genes used': pos_shared_dim1,
            'Number total genes': total_shared_dim1
        },
        {
            'Name': 'Protein abundance',
            'Source': '<a href="https://pax-db.org/download" target="_blank">PaxDb</a>',
            'AUROC': auc_for_paxdb,
            'pvalue': p_for_paxdb,
            'Target genes used': pos_paxdb,
            'Number total genes': total_paxdb
        },
        {
            'Name': 'GC content',
            'Source': '',
            'AUROC': auc_for_gc_content,
            'pvalue': p_for_gc_content,
            'Target genes used': pos_gc_content,
            'Number total genes': total_gc_content
        },
        {
            'Name': 'Inflammatome rank',
            'Source': '<a href="https://www.cell.com/cell-reports/fulltext/S2211-1247(25)01655-9" target="_blank">Cort et al.</a>',
            'AUROC': auc_for_inflammatome,
            'pvalue': p_for_inflammatome,
            'Target genes used': pos_inflammatome,
            'Number total genes': total_inflammatome
        },
    ]
).sort_values('AUROC', ascending=False)
aa_summary_df.index = [""] * len(aa_summary_df)
auc_column = 'AUROC'
auc_format = "{:.2f}"
if single_gene_input:
    aa_summary_df['AUROC'] = aa_summary_df['AUROC'] * 100
    aa_summary_df = aa_summary_df.rename(columns={'AUROC': 'Percentile'})
    auc_column = 'Percentile'
    auc_format = "{:.1f}"
aa_summary_df = aa_summary_df.rename(columns={'pvalue': 'p-value'})
if single_gene_input:
    aa_summary_df = aa_summary_df[
        ['Name', 'Source', 'Target genes used', 'Number total genes', auc_column]
    ]
else:
    aa_summary_df = aa_summary_df[
        ['Name', 'Source', 'Target genes used', 'Number total genes', auc_column, 'p-value']
    ]

aa_summary_display = aa_summary_df.copy()
aa_summary_display['Target genes used'] = aa_summary_display['Target genes used'].astype(int)
aa_summary_display['Number total genes'] = aa_summary_display['Number total genes'].astype(int)
aa_summary_display[auc_column] = aa_summary_display[auc_column].map(
    lambda value: auc_format.format(value) if pd.notna(value) else ""
)
if not single_gene_input:
    aa_summary_display['p-value'] = aa_summary_display['p-value'].map(
        lambda value: f"{value:.2g}" if pd.notna(value) else ""
    )
st.markdown(
    """
<style>
.gene-prop-table-wrapper { width: 100%; overflow-x: scroll; scrollbar-width: auto; }
.gene-prop-table-wrapper::-webkit-scrollbar { height: 13px; }
.gene-prop-table-wrapper::-webkit-scrollbar-thumb { background-color: rgba(0,0,0,0.3); border-radius: 6px; }
.gene-prop-table-wrapper::-webkit-scrollbar-track { background: rgba(0,0,0,0.08); }
table.gene-prop-table { width: 100%; table-layout: auto; }
table.gene-prop-table th, table.gene-prop-table td { white-space: nowrap; }
table.gene-prop-table th { text-align: left; }
</style>
""",
    unsafe_allow_html=True
)
table_html = aa_summary_display.to_html(escape=False, index=False, classes="gene-prop-table")
st.markdown(
    f"<div class=\"gene-prop-table-wrapper\">{table_html}</div>",
    unsafe_allow_html=True
)

if np.isfinite(homlof_fisher_pvalue):
    st.markdown(
        f"Using the {homlof_background_label}, **{homlof_target_hits} of "
        f"{homlof_target_total}** tested input genes are "
        "[genes with homozygous carriers of putative loss-of-function variants]"
        "(https://www.nature.com/articles/s41586-026-10667-5), compared "
        f"with **{homlof_reference_hits} of {homlof_reference_total}** non-target "
        f"background genes. The observed direction is **{homlof_enrichment_direction}** "
        "(two-sided Fisher's exact p-value = "
        f"**{homlof_fisher_pvalue:.2g}**)."
    )
else:
    st.markdown(
        "The enrichment or depletion p-value for genes with homozygous carriers of "
        "putative loss-of-function variants could not be calculated using the "
        f"{homlof_background_label} because there were too few tested input or "
        "background genes."
    )


###classification 
#should be equal to proportions target - needs checking
best_predicted_genes = []

n_splits = 4

if len(target_genes) >= n_splits*2:
    y = proportions['classification_target']
    
    X_proportions = proportions.drop(['classification_target', 'gene_symbol'], axis = 1)
    X_locations = gene_locations.drop(['classification_target', 'gene_symbol'], axis=1)
    y_rankings = table_rankings_intersection['classification_target']
    X_rankings = table_rankings_intersection.drop(['classification_target', 'gene_symbol'], axis=1)

    # Add StandardScaler for gene locations
    scaler = StandardScaler()
    X_locations_scaled = pd.DataFrame(
        scaler.fit_transform(X_locations),
        columns=X_locations.columns,
        index=X_locations.index
    )
    X_rankings_scaled = pd.DataFrame()
    if len(X_rankings) > 0:
        scaler_rankings = StandardScaler()
        X_rankings_scaled = pd.DataFrame(
            scaler_rankings.fit_transform(X_rankings),
            columns=X_rankings.columns,
            index=X_rankings.index
        )


    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=1)
    
    auc_scores = []
    auprc_scores = []
    auc_scores_proportions = []
    auc_scores_locations = []
    auc_scores_rankings = []
    
    with st.spinner('Please wait...'):
      for i, (train_idx, test_idx) in enumerate(skf.split(X_proportions, y)):

          #st.write("fold:" + str(i))
          X_proportions_train = X_proportions.iloc[train_idx, :]
          X_proportions_test = X_proportions.iloc[test_idx, :]
          
          #Gene locations
          X_loc_train = X_locations_scaled.iloc[train_idx, :]
          X_loc_test = X_locations_scaled.iloc[test_idx, :]
          
          y_train = y.iloc[train_idx]
          y_test = y.iloc[test_idx]

          #st.write("fold after mem:" + str(i))
          
          model = LogisticRegression()
          #st.write("fold after init:" + str(i))
  
          gc.collect() 
  
          model.fit(X_proportions_train, y_train)

          # Extract predictions from fitted model
          # probs for classes ordered in same manner as model.classes_
          probas = pd.DataFrame(model.predict_proba(X_proportions_test), columns=model.classes_)

          # Get metrics for proportions + length
          auc = roc_auc_score(y_test, probas[True])
          auc_scores_proportions.append(auc)
          #track the gene with the max prediction for the true class
          best_predicted_genes.append(proportions.iloc[test_idx[probas.idxmax()[True]],:]['gene_symbol'])

          #area under precision recall curve
          auprc = average_precision_score(y_test, probas[True])
          auprc_scores.append(auprc)

          #run the model again with gene locations instead of proportions
          model.fit(X_loc_train, y_train)
          probas = pd.DataFrame(model.predict_proba(X_loc_test), columns=model.classes_)
          auc = roc_auc_score(y_test, probas[True])
          auc_scores_locations.append(auc)

    if len(X_rankings_scaled) > 0 and y_rankings.sum() >= n_splits*2 and (~y_rankings).sum() >= n_splits*2:
      skf_rankings = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=1)
      with st.spinner('Please wait...'):
        for train_idx, test_idx in skf_rankings.split(X_rankings_scaled, y_rankings):
            X_rankings_train = X_rankings_scaled.iloc[train_idx, :]
            X_rankings_test = X_rankings_scaled.iloc[test_idx, :]
            y_rankings_train = y_rankings.iloc[train_idx]
            y_rankings_test = y_rankings.iloc[test_idx]

            model = LogisticRegression()
            model.fit(X_rankings_train, y_rankings_train)
            probas = pd.DataFrame(model.predict_proba(X_rankings_test), columns=model.classes_)
            auc = roc_auc_score(y_rankings_test, probas[True])
            auc_scores_rankings.append(auc)
        
    best_predicted_genes = set(best_predicted_genes)
    top_predicted_hits = best_predicted_genes.intersection(target_genes)
    best_predicted_genes = best_predicted_genes.difference(target_genes)
    
    measures = {
                'number_of_input_genes': len(target_genes),
                'number_of_used_genes': len(target_genes_found),
                'number_of_input_background_genes': len(background_genes),
                'number_of_used_background_genes': len(background_genes_found),
                'area_under_precision_recall_curve_proportions': np.mean(auprc_scores), 
                'AUC proportions': np.mean(auc_scores_proportions),
                'AUC proportions standard dev': np.std(auc_scores_proportions),
                'AUC gene locations': np.mean(auc_scores_locations),
                'AUC all background rankings (no shared dim1)': np.mean(auc_scores_rankings) if len(auc_scores_rankings) else float("nan"),
                'number_of_used_genes_all_rankings_intersection': len(table_rankings_intersection),
                'number_of_target_genes_all_rankings_intersection': int(y_rankings.sum()),
                'AUC proportions p_value versus 0.5' : scistats.ttest_1samp(auc_scores_proportions, 0.5).pvalue, #two sided p-value testing the AUC values against 0.5 expecation
                'AUC locations vrs proportions pvalue' : scistats.ttest_rel(auc_scores_locations, auc_scores_proportions).pvalue}

    st.write("""#### Classification results

To test if the input proteins can be discriminated from the background proteins based on their residue proportions (plus protein coding length) and 
genomic location embeddings. These features were used to train and test a logistic regression model (4 fold cross-validation, 
L2 loss, sklearn default parameters) that attempts to classify proteins as belonging to the input set. Given that the input genes are probably fewer than the background genes, we again report the AUC statistic.

""")
    #st.write('There are  ' + str(X['A']) + ' adenine (A)')

    st.markdown(f"Using amino acid composition and length alone, the average AUC is **{measures['AUC proportions']:.2f}**.")

    st.markdown(f"Using genomic locations alone, the average AUC is **{measures['AUC gene locations']:.2f}**.")

    if np.isfinite(measures['AUC all background rankings (no shared dim1)']):
        st.markdown(
            "Using the intersection of all first-table background rankings "
            "(excluding Shared 1D representation), the average AUC is "
            f"**{measures['AUC all background rankings (no shared dim1)']:.2f}**."
        )
    else:
        st.markdown(
            "Using the intersection of all first-table background rankings "
            "(excluding Shared 1D representation), there are too few genes to run "
            "4-fold cross-validation."
        )

    with st.expander("More statistics from the classification tests"):
        st.write(measures)
    
    #input, X, y, n_jobs, meta_df, name
else:
	st.write("#### Classification results")
	st.write("Too few genes to run classification task - skipping")
