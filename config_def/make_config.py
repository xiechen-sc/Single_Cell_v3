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
    pvalue = 'None'
    splitby = 'None'
    dotsplit = 'False'
    selcet = 'False'
    select_col = 'new_celltype'
    select_lst = ['B_cells',['T_cells','NK']]
    ## 后续引入 mysql 数据库



    # 得到config.yaml
    config_out_file = mkdir(config_out,analysis_type)
    f = open(config_out_file, 'w')
    f.write(f"""
input_seurat: {seurat}  # 输入的 seurat 文件
plot: {plot}  # 可视化方法，可选 ridgeplot,vlnplot,dotplot,featureplot,splitby_featureplot,boxplot
groupby: {groupby} # 分组展示条件，可选clusters、group等
# 下方参数可选
reduct: {reduct}  # 映射降维方式
output: {out}  # 输出目录名
cpu : {cpu}  # CPU 使用数
assay: {assay}  # 正常项目都是 RNA 部分会出现 SCT 或其他
genelist: {genelist}  # 输入的 genlist 文件 
pvalue: {pvalue}  # 为vlnplot和boxplot组间比较，添加p值，可使用all:all 例：--pvalue group1:group2+group2:group3
splitby: {splitby}  # 用来分面绘图的组名,可设为group等
dotsplit: {dotsplit}  # 气泡图是否进行分面展示（默认FALSE），如果为TRUE，根据基因列表有多少列进行gap分隔，并对列名进行简写，如下：
selcet: {selcet}  # 是否 subset 总细胞
select_col: {select_col}  # 选择 subset 的列
select_lst: {select_lst}  # 选择 subset 的列的具体内容, 列表代表 将列表内的细胞合在一起后执行分析
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
    custom_ref = "None"
    white_celltypes = "None"
    assay = 'RNA'
    rerun = 'T'
    extraGene = 'None'
    tissue = 'None'
    celltyping = 'False'
    annolevel = 'single'
    delete_special = 'True'
    # 根据数据库进行修改
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息
    if 'tissue' in project_info:
        tissue = project_info['tissue']
    if 'custom_ref' in project_info:
        custom_ref = project_info['custom_ref']

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
extraGene: {extraGene}  # 额外输入的 marker 基因可视化列表 genelist.txt 若为 None 则不进行核外的 marker 可视化
# 下方内容选择性填写！！！
celltyping: {celltyping}  # 默认不再提供 singleR结果！
delete_special: {delete_special}  # 是否在高变基因中去除列表中的线粒体、热休克、核糖体、解离相关、lncRNA、TR_V_gene和血红蛋白基因（默认去除，仅人小鼠有效）
white_celltypes: {white_celltypes} # 手动指定用哪个细胞类型的经典marker出基因可视化结果,需和上面筛选细胞群一一对应 输入形式为列表，参考：['Macrophages', 'T_cell']
tissue: {tissue} # brain(脑)、Intestinal(肠)、lung(肺)、gastric(胃癌)、tumour(肿瘤). 非必须
custom_ref: {custom_ref}  # 自行指定参考基因组
singleR_rds: {singleR_rds}  #  自动注释参考数据集 如果需要手动指定 请使用绝对路径
annolevel: {annolevel}  # singleR 注释水平  single  main
assay: {assay} # 有时候会用 SCT
rerun: {rerun} # 默认重新寻找高边基因进行聚类
run: {analysis_type} # 这个不要修改
    """)
    f.close()

    print(f'config.yaml 文件已生成至 {config_out_file}')
# 亚群分析系列 old
def get_sub_clusters_old(config_out):
    # 默认变量
    analysis_type2 = 'sub_clusters'
    analysis_type = 'sub_clusters_old'
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
    extraGene = 'None'
    tissue = 'None'
    celltyping = 'False'
    annolevel = 'single'
    delete_special = 'True'
    # 根据数据库进行修改
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息
    if 'tissue' in project_info:
        tissue = project_info['tissue']

    config_out_file = mkdir(config_out,analysis_type2)
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
extraGene: {extraGene}  # 额外输入的 marker 基因可视化列表 genelist.txt 若为 None 则不进行核外的 marker 可视化
# 下方内容选择性填写！！！
celltyping: {celltyping}  # 默认不再提供 singleR结果！
delete_special: {delete_special}  # 是否在高变基因中去除列表中的线粒体、热休克、核糖体、解离相关、lncRNA、TR_V_gene和血红蛋白基因（默认去除，仅人小鼠有效）
tissue: {tissue} # brain(脑)、Intestinal(肠)、lung(肺)、gastric(胃癌)、tumour(肿瘤). 非必须
singleR_rds: {singleR_rds}  #  自动注释参考数据集 如果需要手动指定 请使用绝对路径
annolevel: {annolevel}  # singleR 注释水平  single  main
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
    updata_bynewcelltype = 'True'
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
    top = 20
    analysis_type = 'diff'
    volcano_plot = 'False'
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
volcano_plot: {volcano_plot}  # 是否绘制火山图 默认不出图
top: {top}  # top 绘制热图基因数
run: {analysis_type}   # 这个不要改

""")
    
# reference celltype
def get_singleR(config_out):
    analysis_type = 'singleR'
    seurat = "['seurat.h5seurat']"
    result_perfix = "[]"
    output = 'reference_celltype'
    assay = 'RNA'
    singleR_rds = 'default'
    reduct2 = 'umap'
    species = 'mouse'
    annolevel = 'single'
    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""
seurat: {seurat}  # 输入的 seurat 文件 可一次输入多种
result_perfix: {result_perfix}  # 如果输入多个 seurat 需要每种seurat 对应一种 prefix 防止文件混乱， 如果只做一个亚型注释 可不填此项 ['A','B']
output: {output} # 输出结果目录
assay: {assay} # RNA SCT
singleR_rds: {singleR_rds} # 如果有指定rds 则填写 否则将使用 默认的 rds
reduct2: {reduct2}  # tsne umao
species: {species} # 物种信息
annolevel: {annolevel} # single main
run: {analysis_type} # 这个不要改！
""")
        

def get_enrichment(config_out):
    analysis_type = 'enrichment'
    species = 'mouse'
    input = "['group_CLP-vs-CON-diff-pvalue-0.05-FC-1.5.xls']"
    prefix_lst = '[]'
    analysis_model = "[0]"
    top_n = '50' 
    sort_by = 'gene_diff'
    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息 

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""
# 0:指定基因列表或列表名通配富集;  top200_markers_for_cluster*.xls
# 1:针对流程产出marker_anno基因列表进行clusters拆分并富集;   all_markers_for_each_cluster_anno.xls
# 2:指定差异基因列表富集;  group_CLP-vs-CON-diff-pvalue-0.05-FC-1.5.xls
# 3:指定差异基因目录富集;  /gpfs/oe-scrna/jhyu/project/st/DZOE2023061713-b2-qiankejian-daiwei-cyffpe-m/result/diffexp/
# 4:针对monocle产出module_anno基因列表进行module拆分并富集  pseudotime_heatmap_gene_module_anno.xls
analysis_model: {analysis_model}  # 见上4行注释 必须与input 的 列表元素长度保持一致！ 或者只填写 一种 model
top_n: {top_n}  # 当分析模式为 1 时生效 提供 差异目录结果下的 该文件 SPOCD1_high-vs-low-all_diffexp_genes_anno.xls 每个clusters 取top n个marker基因
sort_by: {sort_by}  # 当分析模式为 1 时生效 提供的 input 文件需要基于哪一列排序 默认 'p-value' 可选 'avg_log2FC' 'q-value' 'gene_diff'
input: {input}  # 输入的 input 文件 可一次输入多种 input  输入可见上方四行注释 可一次填写多个
species: {species} # 物种信息
prefix_lst: {prefix_lst}  # 生成目录前缀，若不填写 则按照input 顺序 写入 0 1 2 3 4......
run: {analysis_type} # 这个不要改！
""")   
        

def get_scenic(config_out):
    analysis_type = 'scenic'
    step1_run = 'False'
    step1_input = 'rds/data_ob_v3.rds'
    species = 'mouse'
    coexMethod = 'top10perTarget'
    step1_outdir = './'
    step2_run = 'False'
    step2_input = 'rds/data_ob_v3.rds'
    result_dir = 'scenic_step1/int/3.4_regulonAUC.Rds'
    rss_rank_top_gene = '3'
    groupby = 'new_celltype'
    sub_seurat = 'False'
    seurat_sub_col = 'new_celltype'
    seurat_sub_col_value = '[]'
    use_color_anno = 'TRUE'
    color_file = ''
    palette = ''
    cluster_n = '4'
    rst_all = '暂无'
    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息 
    if 'scenic' in project_info:
        if 'step1_rds' in project_info['scenic']:
            result_dir = project_info['scenic']['step1_rst'][0]
            rst_all = ','.join(project_info['scenic']['step1_rst'])
    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""
step1_run: {step1_run}  # 是否执行第一步 生成 AUC 活性矩阵
step1_input: {step1_input}  # 输入的 seurat 对象 可以是 rds 也可以是 h5seurat (自动转化)
species: {species}  # 只能是human或者mouse (各种参考基因组版本均可)
coexMethod: {coexMethod}  # 计算调控子共表达的方法（可选：w0.001,w0.005,top50,top50perTarget,top10perTarget,top5perTarget）
step1_outdir: {step1_outdir}  # 活性矩阵输出目录  这个建议别改 改了对脚本也无效
step2_run:  {step2_run}  # 是否执行 RAS CSI 
step2_input:  {step2_input}  # step2 生成 ras csi 结果时 输入的 seurat 对象 可以是 rds 也可以是 h5seurat (自动转化)
result_dir: {result_dir}  # step1 生成的结果文件目录
rss_rank_top_gene: {rss_rank_top_gene}  # RSS rank 绘图结果中加标记的top基因数
groupby: {groupby}  # 图片展示的组别信息
sub_seurat: {sub_seurat}  # 是否对 seurat 对象取子集
seurat_sub_col: {seurat_sub_col}  # 取子集的列名(metadata colnames)
seurat_sub_col_value: {seurat_sub_col_value}  # 取子集的名称，列表内的列表为一次性取出多个， 列表内的字符串元素为分别取子集
# 下方内容选择性填写 建议默认 
use_color_anno: {use_color_anno}  # 是否采用rds中注释的颜色信息,默认为"TRUE"
color_file: {color_file}  # 输入以tab分隔的文件，第一列的列名为metadata列名，第一列为该列元素，第二列为对应的颜色
palette: {palette}  # Get_colors.R 中的离散型色板名 默认"customecol2"
cluster_n:  {cluster_n}  # CSI 聚类数目，第一次跑先设置4，出图后人工判断合适的聚类数，然后重新设置该参数再次运行命令
run: {analysis_type}  # 这个不要改
# 所有 step1 结果可参考下方内容
# {rst_all}
""")

def get_decontX(config_out):
    analysis_type = 'decontX'
    input_seurat = 'seurat.h5seurat'
    threshold = 'NULL'
    reduct = 'NULL'
    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)

    with open(config_out_file,'w')as f:
        f.write(f"""
input_seurat: {input_seurat}  # 输入的 seurat 对象文件 seurat 和 rds 都可以！
threshold: {threshold}  # 污染率阈值 默认不删除，填写范围 0 - 1.0 ，将会删除高于此数值的细胞 
reduct: {reduct}   # 污染可视化，可选用 pca tsne umap ，必须是执行过上述降维后才能使用！   
run: {analysis_type}  # 这个不要改
""")

def get_monocle2(config_out):
    analysis_type = 'monocle2'
    monocle_run = 'FALSE'
    input_seurat = 'seurat.h5seurat'
    assay = 'RNA'
    col_name = 'clusters'
    step1_groupby = 'clusters,sampleid,group'
    sub_seurat = 'False'
    sub_col = 'clusters'
    sub_lst = '1,2,3,4,5'
    output_dir = './monocle2'
    result_rds = 'NULL'
    monocle_rds = 'NULL'
    resolution = '0.4'
    downsample = '30000'
    use_color_anno = 'TRUE'
    color_file = ''
    palette = ''
    cores_use = '8'
    pointsize = '1'
    step2_downstream_run = 'FALSE'
    genelist = 'ordering'
    species = 'human'
    root_state = '1'
    show_branch = 'FALSE'
    branch = 'NULL'
    vis_methods = 'all'
    step2_groupby = 'clusters'
    module_expressplot = 'NULL'
    module_enrichment = 'FALSE'

    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息 

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""
step1_monocle_run: {monocle_run} # 是否执行第一步基础分析
input_seurat: {input_seurat} #输入的 seurat 对象 可以是 rds 也可以是 h5seurat (自动转化)
col_name: {col_name} # 以哪一列筛选高变基因
step1_groupby: {step1_groupby} # 用于作图的分组变量,如clusters,sampleid,group,new_celltype（“State”会默认出图,无需在此填写）
sub_seurat: {sub_seurat} # 是否对输入的seurat对象筛选子集
sub_col: {sub_col} # seurat对象metadata里的分组列名,用于截取seurat的一部分进行分析,如celltype 、clusters、sampleid等
sub_lst: {sub_lst} # 需要截取的分组变量,如有多个可用,分隔
result_rds: {result_rds} # 输入已生成的pseudotime_results.rds, 重新绘图

step2_downstream_run: {step2_downstream_run} # 是否执行monocle下游分析
monocle_rds: {monocle_rds} # 输入已生成的pseudotime_results.rds
species: {species} # 只能是human或者mouse
genelist: {genelist} # 输入基因列表（带表头,表头可随意命名）;或差异基因表格 Diffexp/-vs-.xls;或 ordering:直接对ordering gene作图
root_state: {root_state} # 指定拟时间轨迹起点
show_branch: {show_branch} # 拟时间轨迹图是否展示分支节点。若需要绘制分支热图,可通过此参数查看分支节点。
branch: {branch} # 指定分支节点branchpoint
vis_methods: {vis_methods} # 作图展示形式:heatmap, expressplot, trajectoryplot, treeplot,module, expressplot_line,ridgeplot ,bin, all为全部展示形式都做,module 是按照 module 绘制动力学趋势图;expressplot_line 为基因的分组别动力学趋势图;
step2_groupby: {step2_groupby} # [expressplot、treeplot、expressplot_line、ridgeplot 参数] 作图的分组变量,默认为clusters;当单独绘制山峦图时,可以用逗号分隔多种分组。
module_expressplot: {module_expressplot} # 输入pseudotime_heatmap_gene_module_anno.xls,出具module expressplot,无需和heatmap绑定
module_enrichment: {module_enrichment} #是否进行module基因富集分析
# 下方内容选择性填写 建议默认 
output_dir: {output_dir} # 结果输出目录
assay: {assay} # RNA OR SCT
resolution: {resolution} # 
downsample: {downsample} # 降采样,默认降至30000细胞
cores_use: {cores_use} # 线程数
use_color_anno: {use_color_anno}  # 是否采用rds中注释的颜色信息,默认为"TRUE"
color_file: {color_file}  # 输入以tab分隔的文件,第一列的列名为metadata列名,第一列为该列元素,第二列为对应的颜色
palette: {palette}  # Get_colors.R 中的离散型色板名 默认"customecol2"
run: {analysis_type}  # 这个不要改
        """)

def get_addmodulescore(config_out):
    input  = 'seurat.h5seurat'
    output  = "addmodulescore"
    reduct  = 'umap'
    assay  = "RNA"
    dataslot  = 'counts,data,scale.data' 
    species = 'mouse'
    genelist = 'genelist.txt'
    groupby = 'clusters'
    splitby = 'None'
    fsplitby = 'None'
    pvalue = 'None'
    strict = 'F'
    analysis_type = 'addmodulescore'
    pointsize = '0.5'
    show_box = 'TRUE'
    sub = None
    select = []
    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息 

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""input: {input}  # 输入 seurat 对象文件
species: {species}  # 物种
genelist: {genelist}  # 输入genelist 列表
groupby: {groupby}  # 分组方式
splitby: {splitby}  # 拆分展示小提琴图
fsplitby: {fsplitby}  # 拆分展示umap图
pvalue: {pvalue}  # 添加显著性 all:all 
# 下方内容选择性填写 建议默认 
sub: {sub} # 根据 meta.data 某一列取子集 None 为不取
select: {select}  # 取哪些元素 ['a','b'为] 分别取 a b;[['a','b']] 为一次性取出 a b
reduct: {reduct} # 降维方式 
output: {output}  # 输出目录
show_box: {show_box}  # 小提琴图中是否添加箱图，默认为TRUE
assay: {reduct}  # RNA  SCT
dataslot: {dataslot}  # counts,data,scale.data
strict: {strict}  # 是否使用严格模式筛选gene，默认FALSE
assay: {assay}
pointsize: {pointsize}  # 点的大小 
run: {analysis_type} # 这个不要改！
""")   

#cellchat
def get_cellchat(config_out):
    version = 'cellchat V2'
    input  = 'seurat.h5seurat'
    species = 'mouse'
    output  = "cellchat"
    column4cell = 'clusters'
    groupby = 'None'
    contrast = 'None'
    subsetby = 'None'
    which_cells = 'None'
    topn = '5'
    palette = 'col50'
    rds = 'None'
    strict = 'F'
    cellchatdb = 'None'
    analysis_type = 'cellchat'
    # 数据库交互
    project_info = database_retrieval(config_path=config_out)
    if 'species' in project_info :
        species = project_info['species'] # 更新物种信息 

    config_out_file = mkdir(config_out=config_out,analysis_type=analysis_type)
    with open(config_out_file,'w')as f:
        f.write(f"""
version: {version}  # cellchat 版本 默认V2,可选V1
input: {input}  # 输入 seurat 对象文件(rds文件或h5文件)
species: {species}
output: {output}  # 输出目录
column4cell: {column4cell}  # 用于进行通讯分析的metadata列名
groupby: {groupby}  # 用于比较的分组列名,默认None不分组,如需分组可填group、sampleid等。
contrast: {contrast}  # 默认None不分组,分组则填写case:control;多组用加号连接,如Case1:Control+Case2:Case1
subsetby: {subsetby}  # 需要截取的列名,默认None不截取, (可截取seurat对象,或者cellchat结果对象)
which_cells: {which_cells}  # 需要截取的具体变量(可截取seurat对象,或者cellchat结果对象),如Epicardial_1,Epicardial_2,Epicardial_3
rds: {rds}  # 读取结果rds,重新生成结果,请传入cellchat_list.rds文件
topn: {topn}  # 绘制气泡图筛选的top受配体对数,默认为top5,如果指定了-x默认无筛选,如需筛选可以进行指定。 例：--topn 5
strict: {strict}  # 是否为严格模式，严格模式下，输入的-x参数中文件source,target两列展示严格一对一的受配体对,否则对这两列组合排列生成相互交叉受配体对
palette: {palette}  # 配色
cellchatdb: {cellchatdb}  # cellchat 数据库,默认None
run: {analysis_type} # 这个不要改！
""")