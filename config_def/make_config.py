puple = '\033[35m'
reset = '\033[0m'
cyan = '\033[36m'
yellow = '\033[33m'
# 创建目录
def mkdir(config_out,analysis_type):
    import os
    directory = config_out + '/' + analysis_type
    # 判断目录是否存在
    if not os.path.exists(directory):
        # 创建目录
        os.makedirs(directory)
        print(f"Directory {puple}'{directory}'{reset} created.")
    else:
        print(f"Directory {puple}'{directory}'{reset} already exists.")
    return(directory)

# featureplot 系列绘图
def get_featureplot(config_out):
    
    # 默认变量
    analysis_type = 'featureplot'  # 用于创建目录传入目录名
    seurat = 'seurat.h5seurat'
    out = './featureplot_vlnplot'
    cpu = 2
    assay = 'RNA'
    genelist = 'genelist.txt'
    groupby = 'clusters'
    reduct = 'umap'
    plot = 'vlnplot,featureplot'
    ## 后续引入 mysql 数据库



    # 得到config.yaml
    config_out_file = mkdir(config_out,analysis_type) + '/config.yaml'
    f = open(config_out_file, 'w')
    f.write(f"""run: featureplot  # 这个不要改！
input_seurat: {seurat}  # 输入的 seurat 文件
plot: {plot}  # 可视化方法，可选ridgeplot,vlnplot,dotplot,featureplot,boxplot
groupby: {groupby} # 分组展示条件，可选clusters、group等
# 下方参数可选
reduct: {reduct}  # 映射降维方式
output: {out}  # 输出目录名
cpu : {cpu}  # CPU 使用数
assay: {assay}  # 正常项目都是 RNA 部分会出现 SCT 或其他
genelist: {genelist}  # 输入的 genlist 文件 
    
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')

# 亚群分析系列
def get_sub_clusters():
    pass