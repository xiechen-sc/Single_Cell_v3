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
    init_file = directory + '/.' + analysis_type + "__"
    with open(init_file, "w") as file:
        pass
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
    f.write(f"""
input_seurat: {seurat}  # 输入的 seurat 文件
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
def get_sub_clusters(config_out):
    # 默认变量
    analysis_type = 'sub_clusters'
    seurat = 'seurat.h5seurat'
    species = 'mouse'
    reduct1 = 'pca'
    reduct2 = 'umap'
    batchid = 'batchid'
    resolution = '0.4'
    col_name = 'new_celltype'
    cells = ["T_cells","all",["T_cells","NK"]]
    singleR_rds = 'default'
    assay = 'RNA'
    rerun = 'T'
    # 根据数据库进行修改
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息

    config_out_file = mkdir(config_out,analysis_type)
    f = open(config_out_file, 'w')
    f.write(f"""
seurat: {seurat}  # 输入用于降维的 seurat
species: {species} # 物种
reduct1: {reduct1}  # mnn harmony pca
reduct2: {reduct2} # umap tsne
batchid: {batchid}  # 去批次采用哪一列  部分老师要求使用 sampleid
resolution: {resolution} # 分辨率 T 细胞设置为 0.6, 0.8
col_name: {col_name}  # 根据哪一列选择细胞做亚群分析 如果填写 all 则用所有细胞重新做亚群分析 具体参考下一行
cells: {cells}  # 哪些细胞类型需要做降维 如果需要将两种细胞放在一起降维 可以写成 [["T_cells","NK"],"B_cells"],这样表示将 "T_cells","NK" 两个合在一起降维 ，并对B细胞单独降维
# 下方内容选择性填写！！！
singleR_rds: {singleR_rds}  #  自动注释参考数据集 如果需要手动指定 请使用绝对路径
assay: {assay} # 有时候会用 SCT
rerun: {rerun} # 默认重新寻找高边基因进行聚类
run: {analysis_type} # 这个不要修改
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')

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

# 差异分析 富集分析
def get_diff_enrich(config_out):
    # 默认变量
    seurat = "seurat.h5seurat"
    cell_types = '["1",\'2\']'
    sub_type = "clusters"
    species = 'mouse'
    treat = '["After","a"]'
    control = '["Before" ,"b"]'
    fc = 1.5
    p = 0.05
    vs_type = 'group'
    top = 10
    analysis_type = 'diff'
    # 根据数据库进行修改
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type) 

    # 写入
    with open(config_out_file,'w')as f:
        f.write(f"""
seurat: {seurat} # 输入的 h5seurat 文件
cell_types: {cell_types}  # 哪些clusters 需要做差异分析 列表中每个元素生成一个单独的脚本  如果是 ['all'] 则不提取细胞子集
sub_type: {sub_type}  #  对于上方参数 从metadata中哪一列选择上方列表中的内容
treat: {treat}  # 实验组 组名
control: {control} # 对照组组名  上下一一对应
fc: {fc}  # 差异大小 foldchange
p: {p}  # pvalue 显著性
vs_type: {vs_type}  # 对应上方的 treat control 决定了基于metadata中哪一列选择实验组与对照组
species: {species} # 填写物种
top: {top}  # top 绘制热图基因数
run: {analysis_type}   # 这个不要改

""")
    
