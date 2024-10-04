#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python

def get_top_marker_file(top_n,input,outdir,prefix,sort_by):
    import pandas as pd
    df = pd.read_csv(input,sep='\t')
    df_sorted = df.reindex(df[sort_by].abs().sort_values().index)  # 根据 绝对值大小进行排序
    df_top_n = df_sorted.groupby(sort_by).head(top_n)
    return df_top_n

top_n = 100
input = '/gpfs/oe-scrna/further_analysis/scRNA/Mobi/DZOE2023121146/20240929/3.diff/all-Diffexp/high-vs-low/SPOCD1_high-vs-low-all_diffexp_genes_anno.xls'
outdir = './test'
prefix = 'test_20241004'
sort_by = 'log2FoldChange'

df = get_top_marker_file(top_n,input,outdir,prefix,sort_by)