def get_featureplot(config_out):
    
    # 默认变量
    seurat = 'seurat.h5seurat'
    out = './featureplot_vlnplot'
    cpu = 10
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

    f.write(f"""parameters:
    input_seurat: {seurat}
    output: {out}
    cpu : {cpu}
    assay: {assay}
    genelist: {genelist}
    groupby: {groupby}
    reduct: {reduct}
    plot: {plot}
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')
