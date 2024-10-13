from .base_class import BaseClass
from single_cell_auto.util import *
from single_cell_auto.cmd_module import *
class Sub_Clusters(BaseClass):
    analysis_module = 'sub_clusters'

    def whitelist(self,cell_name,species,tissue,seurat_sub,reduct2,cell_name_out):  # 白名单相关 
        database = '/data/database/sc_subtype_refmarker/'
        whitelist_file = database + 'whitelist.xls'
        if cell_name == 'all':
            print("既未指定细胞类型，又无法获取有效的细胞类型信息，不执行基因可视化。")
            return 'None'
        elif species not in Frequent_species():
            print("非常规物种不适用默认基因列表")
            return 'None'
        elif tissue not in  ['brain','Intestinal','lung','gastric','tumour']:
            print("未提供有效组织类型信息，提供非组织类型marker基因可视化，若无请核查。")
        if type(cell_name) == list:
            all_celltype = [i for i in cell_name]
            cell_name_out = "_".join(cell_name)
        else:
            all_celltype = [cell_name]
            cell_name_out = cell_name
        cell_name_out = cell_name_normalization(cell_name_out)
        f = open(whitelist_file,'r')
        lines = f.readlines()
        f.close()
        # 将白名单 转化为字典 如果匹配到值 则返回键！
        col_num = len(lines[0].replace('\n','').split('\t'))
        info = {i:lines[0].strip().split('\t')[i] for i in range(col_num)}
        line_counts = 0
        celltype = {}  
        for line in lines:  # 将白名单 转化为字典 如果匹配到值 则返回键！
            line = line.replace('\n','').split('\t')   
            if line_counts == 0:
                line_counts = 1
                for j in line:
                    celltype[j] = []
            else:
                for k in range(col_num):
                    celltype[info[k]].append(line[k])

        use_celltype = []  # celltype 包含的所有白名单细胞
        
        for j in all_celltype:
            for k in celltype:
                if j in celltype[k]:
                    use_celltype.append(k)
        
        use_celltype = [(database + i) for i in use_celltype] # 将所有匹配到的目录保存

        import os
        use_celltype_list = []
        for directory_path in use_celltype:
            marker_file = os.listdir(directory_path)
            marker_file = [os.path.join(directory_path, d) for d in marker_file]
            for f in marker_file:
                use_celltype_list.append(f)
        

        import re
        if tissue == 'None':
            tissue_pattern = r"reference_marker_\d+_[a-zA-Z]+\.xls$"
            cmd1 = ""
            for mkfile in use_celltype_list:
                if not bool(re.findall(tissue_pattern,mkfile)):
                    outcellname = "_".join(cell_name)
                    cmd1 += f"""
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
  -i {seurat_sub}  \\
  -f h5seurat \\
  -o sub_{cell_name_out}/featureplot_vlnplot \\
  -j 10 \\
  --assay RNA \\
  --dataslot data \\
  visualize \\
  -x {mkfile} \\
  -g clusters \\
  --reduct {reduct2} \\
  -m vlnplot,featureplot,dotplot \\
  --vcolors customecol2 \\
  --ccolors spectral \\
  --pointsize 0 \\
  --dodge F                    
"""
                    
            return cmd1
        else:
            cmd1 = ""
            for mkfile in use_celltype_list:
                if bool(re.findall(tissue,mkfile)):
                    cmd1 += f"""
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
  -i {seurat_sub}  \\
  -f h5seurat \\
  -o sub_{cell_name_out}/featureplot_vlnplot \\
  -j 10 \\
  --assay RNA \\
  --dataslot data \\
  visualize \\
  -x {mkfile} \\
  -g clusters \\
  --reduct {reduct2} \\
  -m vlnplot,featureplot,dotplot \\
  --vcolors customecol2 \\
  --ccolors spectral \\
  --pointsize 0 \\
  --dodge F                    
"""
            return cmd1
        
    def get_cell_type(cells):
        pass
    def get_script(self):
        annolevel = self.annolevel
        celltyping = self.celltyping
        tissue = self.tissue
        extraGene = self.extraGene
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
        delete_special = self.delete_special

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
            anno_dir = anno
            anno = anno + 'annotation/gene_annotation.xls'

        
            
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        
        
        # 处理 cell
        for j in cells:
            if type(j) == list:
                str_list = [str(k) for k in j]
                cell_name = "_".join(str_list)
                cell_type = ",".join(["\\'" + k + "\\'" for k in str_list])
                cell_name_raw = str_list
            elif type(j) == str:
                cell_name_raw = str(j)
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
--prefix {prefix} \\
"""
            if cell_name != 'all':
                cmd += f'--predicate  "{col_name} %in% c({cell_type})"   \\\n'
            cmd += f"bclust   \\\n--reduct1 {reduct1}  \\"
            if reduct1 == 'mnn':
                cmd += f"\n--batchid {batchid} \\\n--components 10  \\"
            elif reduct1 == "pca,harmony":
                cmd += f"\n--batchid {batchid} \\\n-t 20 \\\n-y 30 \\"
            elif reduct1 == 'pca':
                pass
            else:
                jinggao("reduct type error,only pca mnn harmony!")
                exit()
            cmd += f"""
--reduct2 {reduct2}   \\
--clusteringuse snn  \\
--resolution {resolution}   \\
--rerun {rerun}   \\
--pointsize  0.5  \\
--palette customecol2"""     
            if delete_special:  # 判断是否在人和小鼠中
                cmd += f""" \\
--ref {anno_dir}

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
            if bool(celltyping):
                if singleR_rds == 'default':
                    try:
                        singleR_rds = species_info['singleR']['default']
                    except KeyError:
                        singleR_rds = '# 请手动填写！！！'
                        jinggao(f'{species} 的 singleR参考注释文件 在数据库中不存在 脚本  singleR 部分 请手动填写参考数据集rds 或删除不做！！')
                    except TypeError:
                        singleR_rds = '# 请手动填写！！！'
                        jinggao(f'{species} 在数据库中不存在 脚本  singleR 部分 请手动填写参考数据集rds 或删除不做！！')
                    else:
                        pass
                else:
                    pass
                species_input2singleR = re.sub('_.*','',species)               
                cmd += cmd_singleR(seurat=seurat_sub,
                                   output=f"sub_{cell_name}/Reference_celltype",
                                   assay=assay,
                                   singleR_rds=singleR_rds,
                                   reduct2=reduct2,
                                   species=species_input2singleR,
                                   annolevel=annolevel) 

            # genelist vis marker 
            if extraGene != 'None':
                cmd += f"""
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
  -i {seurat_sub}  \\
  -f h5seurat \\
  -o sub_{cell_name}/featureplot_vlnplot \\
  -j 10 \\
  --assay RNA \\
  --dataslot data \\
  visualize \\
  -x {extraGene} \\
  -g clusters \\
  --reduct {reduct2} \\
  -m vlnplot,featureplot,dotplot \\
  --vcolors customecol2 \\
  --ccolors spectral \\
  --pointsize 0 \\
  --dodge F
"""
            cmd1 = self.whitelist(cell_name=cell_name_raw,species=species,tissue=tissue,seurat_sub=seurat_sub,reduct2=reduct2,cell_name_out=cell_name)
            if cmd1 != 'None':
                cmd += cmd1
                
            with open(out_script,"w") as f:
                f.write(cmd)
            print(f"脚本 {out_script} 已生成") 
        #### 物种信息保存至数据库
        db_update_bg = self.pjif  
        db_update_bg['species'] = species
        self.update_info_bag = db_update_bg
        
            
            
            
if __name__ == '__main__':
    pass