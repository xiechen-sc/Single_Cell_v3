from .base_class import BaseClass
from single_cell_auto.util import *

class Modified_cell_type(BaseClass):
    analysis_module = 'modified_cell_type'

    def get_script(self):
        seurat = self.input
        output = self.output
        updata = self.updata
        updata_bynewcelltype = self.updata_bynewcelltype  
        type_name = self.type_name
        newseurat = self.newseurat
        species = self.species 
        custom_ref = self.custom_ref
        if updata:
            bl = 'TRUE'
        else:
            bl = 'FALSE'
        newcelltype_file = self.Modified_file
        Modified_col = self.Modified_col
        reduct = self.reduct
        newcelltype_file_type = newcelltype_file.split('.')
        newcelltype_file_type = newcelltype_file_type[len(newcelltype_file_type)-1]
        species_info = get_species_info(species=species)
        if custom_ref != "None":
            anno_dir = custom_ref
            anno = anno_dir + '/annotation/gene_annotation.xls'
        else:
            try:
                anno = species_info['anno']
            except KeyError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
            except TypeError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 在数据库中不存在 请手动填写！')
            else:
                anno = anno + 'annotation/gene_annotation.xls'
        if newcelltype_file_type == 'tsv':
            fgf = 'F'
        elif newcelltype_file_type == 'csv':
            fgf = 'T'
        else:
            print('细胞文件名后缀只能是 tsv 或 csv，请检查 config 文件中 Modified_file 参数！')
            exit(1)
        out_script = f'{self.outdir}/cmd_modified_cell_type.sh'
        #### 物种信息保存至数据库
        database_add(config_path=out_script,config_info={'species':species})  
        with open(out_script,"w") as f:
            f.write(f"""set -e
module purge && module load OESingleCell/3.0.d
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat} \\
-f h5seurat \\
-o {output} \\
-d h5seurat \\
--update {bl} \\
--assay RNA \\
--dataslot counts,data,scale.data  \\
changecelltype \\
-c {newcelltype_file} \\
-C {Modified_col} \\
--palette customecol2 \\
--reduct {reduct} \\
-b {fgf}
                """)

            if updata_bynewcelltype:
                seurat=newseurat
                f.write(f"""
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i  {seurat} \\
-f h5seurat \\
-o ./ \\
--assay RNA \\
--dataslot data \\
summarize \\
--reduct umap \\
--palette customecol2 \\
-c {type_name} \\
-b sampleid,group \\
--pointsize 0.5 \\
--dosummary T


# 1.鉴定marker
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat} \\
-f h5seurat \\
-o Marker \\
--assay RNA \\
--dataslot data,counts \\
-j 10 \\
findallmarkers \\
-c 2 \\
-N 10 \\
-k 1 \\
-p 0.05 \\
-s F \\
-e presto \\
-n {type_name}

#2.可视化+anno
# marker热图
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \\
-i {seurat} \\
-f h5seurat \\
-o ./Marker \\
-t 10 \\
--assay RNA \\
--slot data,scale.data \\
heatmap \\
-l Marker/top10_markers_for_each_cluster.xls \\
-c gene_diff \\
-n 10 \\
-g {type_name} \\
--group_colors customecol2 \\
--sample_ratio 0.8 \\
--style seurat

# featureplot小提琴图
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat}  \\
-f h5seurat \\
-o ./Marker \\
-j 10 \
--assay RNA \\
--dataslot data \\
visualize \
-l Marker/top10_markers_for_each_cluster.xls \
-g {type_name} \\
--reduct umap \\
--topn  10  \\
--topby gene_diff \\
-m vlnplot,featureplot \\
--vcolors customecol2 \\
--ccolors spectral \\
--pointsize 0.3 \\
--dodge F

#anno


Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \\
-g Marker/all_markers_for_each_cluster.xls \\
--anno {anno}  # 根据物种修改
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \\
-g Marker/top10_markers_for_each_cluster.xls \\
--anno {anno}  # 根据物种修改
""")

        print(f"脚本 {out_script} 已生成")
        db_update_bg = self.pjif  
        db_update_bg['species'] = species
        db_update_bg['custom_ref'] = custom_ref
        self.update_info_bag = db_update_bg