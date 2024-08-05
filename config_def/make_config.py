from single_cell_auto.util import *
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
    return(directory + '/config.yaml')

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
    config_out_file = mkdir(config_out,analysis_type)
    f = open(config_out_file, 'w')
    f.write(f"""input_seurat: {seurat}  # 输入的 seurat 文件
plot: {plot}  # 可视化方法，可选ridgeplot,vlnplot,dotplot,featureplot,boxplot
groupby: {groupby} # 分组展示条件，可选clusters、group等
# 下方参数可选
reduct: {reduct}  # 映射降维方式
output: {out}  # 输出目录名
cpu : {cpu}  # CPU 使用数
assay: {assay}  # 正常项目都是 RNA 部分会出现 SCT 或其他
genelist: {genelist}  # 输入的 genlist 文件 
run: featureplot  # 这个不要改！
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')

# 亚群分析系列
def get_sub_clusters():
    pass

# 修改细胞类型
def get_modified_cell_type(config_out):
    # 默认变量
    analysis_type = 'modified_cell_type'  # 用于创建目录传入目录名
    seurat = 'seurat.h5seurat'
    updata = 'False'
    output = 'newcelltype'
    Modified_file = 'newcelltype.tsv '
    Modified_col = 'clusters'
    reduct = 'umap'
    updata_bynewcelltype = 'False'
    newseurat = 'newcelltype/seurat.h5seurat'
    type_name = 'new_celltype'
    species = 'mouse'
    # 根据数据库进行修改
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息
        



    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type) 
    with open(config_out_file,'w')as f:
        f.write(f"""
input : {seurat}  #输入的 seurat 文件
updata: {updata}  # 是否覆盖原文件 False 则在输出目录生成注释后的seurat
output: {output}  # 结果输出目录
Modified_file: {Modified_file}  # 用来修改细胞类型的文件
Modified_col: {Modified_col}  #  用来修改的列
reduct: {reduct}  # 降维方法
updata_bynewcelltype: {updata_bynewcelltype}  # 是否更新后续基于 newcelltype的分析 包含 cor vis marker 分析
newseurat: {newseurat}  # 更新细胞类型后 新细胞类型的列
type_name: {type_name}  # 
species: {species}  # 物种 
run: {analysis_type}  # 这个不要改
""")
    print(f'config.yaml 文件已生成至 {config_out_file}')