from .base_class import BaseClass
from single_cell_auto.util import *
class Sub_Clusters(BaseClass):
    analysis_module = 'sub_clusters'

    def whitelist():  # 白名单相关
        pass
    def get_cell_type(cells):
        pass
    def get_script(self):
        seurat = self.seurat
        species = self.species
        reduct1 = self.reduct1
        reduct2 = self.reduct2  
        batchid = self.batchid
        resolution = self.resolution
        col_name = self.col_name 
        cells = self.cells
        assay = self.assay
        rerun = self.rerun
        singleR_rds = self.singleR_rds

        species_info = get_species_info(species=species)
        get_anno = True
        
        try:
            anno = species_info['anno'] 
        except KeyError:
            anno = '# 请手动填写！！！'
            jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
            get_anno = False
        except TypeError:
            anno = '# 请手动填写！！！'
            jinggao(f'{species} 在数据库中不存在 请手动填写！')
            get_anno = False
        else:
            anno = anno + 'annotation/gene_annotation.xls'

        if singleR_rds == 'default':
            try:
                singleR_rds = species_info['singleR']['default']
            except KeyError:
                singleR_rds = '# 请手动填写！！！'
                jinggao(f'{species} 的 singleR参考注释文件 在数据库中不存在 脚本 {out_script} singleR 部分 请手动填写参考数据集rds 或删除不做！！')
            except TypeError:
                singleR_rds = '# 请手动填写！！！'
                jinggao(f'{species} 在数据库中不存在 脚本 {out_script} singleR 部分 请手动填写参考数据集rds 或删除不做！！')
            else:
                anno = anno + 'annotation/gene_annotation.xls'
        else:
            singleR_rds = self.singleR_rds
            
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        #### 物种信息保存至数据库
        database_add(config_path=out_script,config_info={'species':species}) 
        
        # 处理 cell
        for j in cells:
            if type(j) == list:
                str_list = [str(k) for k in j]
                cell_name = "_".join(str_list)
                cell_type = ",".join(["\\'" + k + "\\'" for k in str_list])
            elif type(j) == str:
                cell_name = str(j)
                cell_type = "\\'" + str(j) + "\\'"
            else:
                exit("config.yaml 文件中 sub_clusters cells 填写格式错误 请查看注释信息！")
            cell_name = cell_name_normalization(cell_name)
            out_script = f'{self.outdir}/cmd_sub_{col_name}_{cell_name}.sh'
            prefix = col_name + "_" + cell_name
            seurat_sub=f"sub_{cell_name}/Clustering/{prefix}.h5seurat"
            # 定义空命令
            if reduct1 == 'harmony':
                reduct1 = "pca,harmony"
            cmd = "set -e\nmodule purge && module load OESingleCell/3.0.d\n"
            cmd = cmd + f"""# reduct
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  \\
-i  {seurat}  \\
-f h5seurat  \\
-o sub_{cell_name}/Clustering  \\
-d h5seurat   \\
--assay {assay}  \\
--dataslot counts,data,scale.data   \\
--update F   \\
--prefix {prefix}
"""
            if cell_name != 'all':
                cmd += f'--predicate  "{col_name} %in% c({cell_type})"   \\\n'
            cmd += f"bclust   \\\n--reduct1 {reduct1}  \\\n"
            if reduct1 == 'mnn':
                cmd += f"--batchid {batchid} \\\n--components 10  \\"
            elif reduct1 == "pca,harmony":
                cmd += f"--batchid {batchid} \\\n-t 20 \\\n-y 30 \\"
            else:
                jinggao("reduct type error,only pca mnn harmony!")
                exit()
            cmd += f"""
--reduct2 {reduct2}   \\
--clusteringuse snn  \\
--resolution {resolution}   \\
--rerun {rerun}   \\
--pointsize  0.5  \\
--palette customecol2

"""
            
            # vis by clusters
            cmd += f"""
# vis by clusters
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat_sub} \\
-f h5seurat \\
-o sub_{cell_name}/Clustering \\
--assay {assay} \\
--dataslot data \\
summarize \\
--reduct {reduct2} \\
--palette customecol2 \\
-c clusters \\
-b sampleid,group \\
--pointsize 0.5 \\
--dosummary T
"""
            # cor plot
            cmd += f"""
# cor analysis
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \\
-i {seurat_sub}  \\
-f h5seurat \\
-o sub_{cell_name}/Clustering/clusters_correlation \\
-t 6 \\
--assay {assay} \\
--slot data \\
--reduct {reduct2} \\
coefficient \\
-g clusters
"""
            # marker
            cmd += f"""
# marker
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat_sub}  \\
-f h5seurat \\
-o sub_{cell_name}/Marker \\
--assay {assay} \\
--dataslot data,counts \\
-j 10 \\
findallmarkers \\
-c 2 \\
-N 10 \\
-k 1 \\
-p 0.05 \\
-s F \\
-e presto \\
-n clusters
"""
            # marker vis + anno
            cmd += f"""

# marker vis + anno
# marker热图
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \\
-i {seurat_sub}  \\
-f h5seurat \\
-o sub_{cell_name}/Marker \\
-t 10 \\
--assay {assay} \\
--slot data,scale.data \\
heatmap \\
-l sub_{cell_name}/Marker/top10_markers_for_each_cluster.xls \\
-c gene_diff \\
-n 10 \\
-g clusters \\
--group_colors customecol2 \\
--sample_ratio 0.8 \\
--style seurat

Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
-i {seurat_sub}   \\
-f h5seurat \\
-o sub_{cell_name}/Marker \\
-j 10 \\
--assay {assay} \\
--dataslot data \\
visualize \\
-l sub_{cell_name}/Marker/top10_markers_for_each_cluster.xls \\
-g clusters \\
--reduct {reduct2} \\
--topn  10  \\
--topby gene_diff \\
-m vlnplot,featureplot \\
--vcolors customecol2 \\
--ccolors spectral \\
--pointsize 0.3 \\
--dodge F

Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \\
-g sub_{cell_name}/Marker/all_markers_for_each_cluster.xls \\
--anno {anno}  # 根据物种修改
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \\
-g sub_{cell_name}/Marker/top10_markers_for_each_cluster.xls \\
--anno {anno}  # 根据物种修改
"""
            # singleR anno
            cmd += f"""
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  \\
-i {seurat_sub}  \\
-f h5seurat \\
-o sub_{cell_name}/Reference_celltype \\
-d h5seurat \\
--update T \\
--assay {assay} \\
--dataslot counts \\
celltyping \\
-r {singleR_rds} \\
--annolevel single \\
--usecluster F \\
--demethod classic \\
--pointsize 0.3 \\
-n 25 \\
--reduct {reduct2} \\
--species {species}
"""

            with open(out_script,"w") as f:
                f.write(cmd)
            print(f"脚本 {out_script} 已生成") 