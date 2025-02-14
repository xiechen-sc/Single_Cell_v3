from single_cell_auto.util import jinggao
# singleR 注释
def cmd_singleR(seurat,output,assay,singleR_rds,reduct2,species,annolevel):
    cmd = f"""set -e
module purge && module load OESingleCell/3.0.d
# reference celltype
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  \\
-i {seurat}  \\
-f h5seurat \\
-o {output} \\
-d h5seurat \\
--update T \\
--assay {assay} \\
--dataslot counts \\
celltyping \\
-r {singleR_rds} \\
--annolevel {annolevel} \\
--usecluster F \\
--demethod classic \\
--pointsize 0.3 \\
-n 25 \\
--reduct {reduct2} \\
--species {species}
"""
    return cmd


def volcano(input,pvalue,log2fc,output):
    cmd = f"""
set -e
module purge && module load OESingleCell/3.0.d
Rscript /gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/volcanoplot/volcano.r \\
-i {input} \\
-p {pvalue} \\
-f {log2fc} \\
-o {output}
"""
    return cmd


def cmd_addmodulescore(self):
    
    analysis_module = self.run
    input = self.input # 输入 seurat 对象文件
    species = self.species  # 物种
    genelist = self.genelist  # 输入genelist 列表
    groupby = self.groupby  # 分组方式
    splitby = self.splitby  # 拆分展示小提琴图
    fsplitby = self.fsplitby  # 拆分展示umap图
    pvalue = self.pvalue
  
    # 下方内容选择性填写 建议默认 
    show_box = self.pvalue
    reduct = self.reduct # 降维方式 
    output = self.output  # 输出目录
    assay = self.reduct  # RNA  SCT
    dataslot = self.dataslot  # counts,data,scale.data
    strict = self.strict  # 是否使用严格模式筛选gene，默认FALSE
    assay = self.assay
    pointsize = self.pointsize

    try:
        anno = self.pjif['species']
    except KeyError:
        anno = '# 请手动填写！！！'
        jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
    except TypeError:
        anno = '# 请手动填写！！！'
        jinggao(f'{species} 在数据库中不存在 请手动填写！')

    cmd = f"set -e\nmodule purge && module load OESingleCell/3.0.d\n"
    cmd += f"""Rscript /gpfs/oe-scrna/chenhaoruo/script/new-sctool/sctool \\
    --input {input} \\
    --reduct {reduct} \\
    --output {output} \\
    --assay {assay} \\
    --dataslot {dataslot} \\
    --anno {anno} \\
    sc_addmodulescore  \\
        -x {genelist} \\
        -g {groupby} \\
"""
    if not (splitby == 'None'):
        cmd += f'    	--splitby {splitby} \\\n'
    if not (fsplitby == 'None'):
        cmd += f'    	--fsplitby {fsplitby} \\\n'
    if not (show_box  == 'None'):
        cmd += f'    	--show_box {show_box} \\\n'
    if not (strict == 'None'):
        cmd += f'    	--strict {strict} \\\n'
    if not (pvalue == 'None'):
        cmd += f'    	--pvalue {pvalue} \\\n'
    cmd += f'        --pointsize {pointsize} '
    return cmd