def get_featureplot(config_out):
    
    # 默认变量
    seurat = 'seurat.h5seurat'
    out = './featureplot_vlnplot'
    cpu = 2
    assay = 'RNA'
    genelist = 'genelist.txt'
    groupby = 'clusters'
    reduct = 'umap'
    plot = 'vlnplot,featureplot'
    ## 后续引入 mysql 数据库


    ##
    # 得到config.yaml
    config_out_file = config_out + '/config.yaml'
    f = open(config_out_file, 'w')

    f.write(f"""run: featureplot  # 这个不要改！
input_seurat: {seurat}
plot: {plot}
groupby: {groupby}
reduct: {reduct}
# 下方参数可选
output: {out} 
cpu : {cpu}
assay: {assay}
genelist: {genelist}
    
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')
