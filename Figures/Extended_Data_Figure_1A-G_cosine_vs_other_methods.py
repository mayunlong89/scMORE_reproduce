import matplotlib.pyplot as plt
import IPython
import scanpy as sc
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import cosg
import numpy as np
plt.rcParams.update({'axes.labelsize' : 'large'})
matplotlib.use('TkAgg')
groupby='Seurat_cell_type'
os.chdir(r'D:\生信\ctDRTF\10X_brain')

#Save plot data from Dotplot
def extract_dotplot_data(dotplot_obj, df_tmp=None, output_path="dotplot_data.csv"):

    color_df = dotplot_obj.dot_color_df.reset_index()
    color_long = color_df.melt(
        id_vars='Seurat_cell_type',
        var_name='Gene',
        value_name='Mean_expression_in_group'
    ).rename(columns={'Seurat_cell_type': 'Cell_type'})

    size_df = dotplot_obj.dot_size_df.reset_index()
    size_long = size_df.melt(
        id_vars='Seurat_cell_type',
        var_name='Gene',
        value_name='Fraction_of_cells_in_group(%)'
    ).rename(columns={'Seurat_cell_type': 'Cell_type'})

    combined_df = pd.merge(
        color_long,
        size_long,
        on=['Cell_type', 'Gene'],
        how='inner'
    )

    if df_tmp is not None:
        target_order = df_tmp.index.tolist()
        valid_order = [ct for ct in target_order if ct in combined_df['Cell_type'].unique()]
        combined_df['Cell_type'] = pd.Categorical(
            combined_df['Cell_type'],
            categories=valid_order,
            ordered=True
        )
        combined_df = combined_df.sort_values(by=['Gene', 'Cell_type']).reset_index(drop=True)

    combined_df.to_csv(output_path, index=False)
    print(f"Data write to: {output_path}")

    return combined_df

#Load pbmc 2k data and preprocess
adata = sc.read_h5ad('adata.h5ad')
sc.tl.pca(adata)
sc.pl.pca(adata, color = "Seurat_cell_type")

sc.tl.rank_genes_groups(adata,groupby='Seurat_cell_type',method = 'wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata,groupby='Seurat_cell_type',cmap='Spectral_r',
                                standard_scale='var',n_genes=3)

sc.tl.dendrogram(adata,groupby=groupby,use_rep='X_pca')

#Run Wilcoxon-test
df_tmp=pd.DataFrame(adata.uns['rank_genes_groups']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'Wilcoxon-test',save = 'Wilcoxon-test.pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                            title='Wilcoxon-test', save='Wilcoxon-test.pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='Wilcoxon-test.csv')

#Run Wilcoxon-test(TIE)
sc.tl.rank_genes_groups(adata,groupby='Seurat_cell_type',method = 'wilcoxon',tie_correct= True)
df_tmp=pd.DataFrame(adata.uns['rank_genes_groups']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'Wilcoxon-test(TIE)',save = 'Wilcoxon-test(TIE).pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'Wilcoxon-test(TIE)',save = 'Wilcoxon-test(TIE).pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='Wilcoxon-test(TIE).csv')

#Run T-test
sc.tl.rank_genes_groups(adata,groupby='Seurat_cell_type',method = 't-test')
df_tmp=pd.DataFrame(adata.uns['rank_genes_groups']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'T-test',save = 'T-test.pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'T-test',save = 'T-test.pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='T-test.csv')

#Run T-test Overestim Var
sc.tl.rank_genes_groups(adata,groupby='Seurat_cell_type',method = 't-test_overestim_var')
df_tmp=pd.DataFrame(adata.uns['rank_genes_groups']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'T-test Overestim Var',save = 't-test_overestim_var.pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'T-test Overestim Var',save = 'T-test Overestim Var.pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='T-test Overestim Var.csv')

#Run Mean Gene Expression
groupby = 'Seurat_cell_type'

mean_exp = pd.DataFrame(
    index=adata.var_names,
    columns=adata.obs[groupby].unique()
)

for ct in mean_exp.columns:
    # 提取该细胞类型的所有细胞
    ct_mask = adata.obs[groupby] == ct
    # 计算该细胞类型内基因的平均表达量
    mean_exp[ct] = np.array(adata[ct_mask, :].X.mean(axis=0)).flatten()

# 对每个细胞类型的基因按平均表达量降序排序，并提取 top3
top3_genes = {}
for ct in mean_exp.columns:
    top3_genes[ct] = mean_exp[ct].sort_values(ascending=False).head(3).index.tolist()

marker_genes_list = top3_genes
if f'dendrogram_{groupby}' not in adata.uns:
    sc.tl.dendrogram(adata, groupby=groupby, use_rep='X_pca')

sc.pl.dotplot(
    adata,
    marker_genes_list,
    groupby=groupby,
    dendrogram=True,
    swap_axes=False,
    standard_scale='var',
    cmap='Spectral_r',
    title='Mean Gene Expression',
    save='Mean_gene_expression.pdf'
)

dotplot_obj = sc.pl.dotplot(
    adata,
    marker_genes_list,
    groupby=groupby,
    dendrogram=True,
    swap_axes=False,
    standard_scale='var',
    cmap='Spectral_r',
    title='Mean Gene Expression',
    save='Mean_gene_expression.pdf',
    return_fig=True
)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='Mean_gene_expression.csv')

#Run COSG
cosg.cosg(adata,
    key_added='cosg',
    mu=100,
    expressed_pct=0.1,
    remove_lowly_expressed=True,
     n_genes_user=2000,
               groupby=groupby)

df_tmp=pd.DataFrame(adata.uns['cosg']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}

sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'COSG',save = 'COSG.pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'COSG',save = 'COSG.pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='COSG.csv')

#Run Logistic Regression
sc.tl.rank_genes_groups(adata,groupby='Seurat_cell_type',method = 'logreg',max_iter=10000)
df_tmp=pd.DataFrame(adata.uns['rank_genes_groups']['names'][:3,]).T
df_tmp=df_tmp.reindex(adata.uns['dendrogram_'+groupby]['categories_ordered'])
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'Logistic regression',save = 'Logistic_egression.pdf')
dotplot_obj = sc.pl.dotplot(adata, marker_genes_list,groupby=groupby,dendrogram=True,swap_axes=False,standard_scale='var',cmap='Spectral_r',
                   title = 'Logistic regression',save = 'Logistic_egression.pdf',return_fig=True)
extract_dotplot_data(dotplot_obj,df_tmp,output_path='Logistic_egression.csv')
