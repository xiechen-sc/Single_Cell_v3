from .base_class import BaseClass
from single_cell_auto.util import *
from single_cell_auto.cmd_module import *
class Sub_Clusters(BaseClass):
    analysis_module = 'sub_clusters'

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
        white_celltypes = self.white_celltypes
        assay = self.assay
        rerun = self.rerun
        singleR_rds = self.singleR_rds
        custom_ref = self.custom_ref
        delete_special = self.delete_special

        species_info = get_species_info(species=species)
        get_anno = True
        
        if custom_ref != "None":
            anno_dir = custom_ref
            anno = anno_dir + '/annotation/gene_annotation.xls'
        else:
            try:
                anno = species_info['anno'] 
            except KeyError:
                anno = '# 请手动填写！！！'
                anno_dir = anno
                jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
                get_anno = False
            except TypeError:
                anno = '# 请手动填写！！！'
                anno_dir = anno
                jinggao(f'{species} 在数据库中不存在 请手动填写！')
                get_anno = False
            else:
                anno_dir = anno
                anno = anno + 'annotation/gene_annotation.xls'

        
            
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        
        
        # 处理 cell
        num = 0
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
            cmd += self.add_cmd_row(f'Rscript /gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/Subtype_new_sctool.R')
            cmd += self.add_cmd_row(f'-i {seurat}')
            cmd += self.add_cmd_row(f'--assay {assay}')
            if cell_name != 'all':
                cmd += self.add_cmd_row(f'--predicate  "{col_name} %in% c({cell_type})" ')
            cmd += self.add_cmd_row(f'--reduct1 {reduct1}')
            if reduct1 == 'mnn':
                cmd += self.add_cmd_row(f'--batchid {batchid}')
            elif reduct1 == "pca,harmony":
                cmd += self.add_cmd_row(f'--batchid {batchid}')
            elif reduct1 == 'pca':
                pass
            else:
                jinggao("reduct type error,only pca mnn harmony!")
                exit()
            cmd += self.add_cmd_row(f'--reduct2 {reduct2}')
            cmd += self.add_cmd_row(f'--resolution {resolution}') 
            if bool(celltyping):
                cmd += self.add_cmd_row(f'-m clustering,marker,celltyping,marker_vis')
            else:
                cmd += self.add_cmd_row(f'-m clustering,marker,marker_vis')                
            species_input2singleR = re.sub('_.*','',species)  
            cmd += self.add_cmd_row(f'-s {species_input2singleR}')
            cmd += self.add_cmd_row(f'--anno {anno}')

            # singleR anno
            if bool(celltyping):
                if singleR_rds == 'default':
                    try:
                        singleR_rds = species_info['singleR']['default']
                        cmd += self.add_cmd_row(f'-r {singleR_rds}')             
                        cmd += self.add_cmd_row(f'--annolevel {annolevel}')   
                    except KeyError:
                        singleR_rds = '# 请手动填写！！！'
                        jinggao(f'{species} 的 singleR参考注释文件 在数据库中不存在 脚本  singleR 部分 请手动填写参考数据集rds 或删除不做！！')
                    except TypeError:
                        singleR_rds = '# 请手动填写！！！'
                        jinggao(f'{species} 在数据库中不存在 脚本  singleR 部分 请手动填写参考数据集rds 或删除不做！！')
                    else:
                        pass
                else:
                    cmd += self.add_cmd_row(f'-r {singleR_rds}')             
                    cmd += self.add_cmd_row(f'--annolevel {annolevel}')         

            # genelist vis marker 
            if white_celltypes != 'None':
                if type(white_celltypes) == list:
                    if len(white_celltypes[num]) != 0:
                        white_celltype = white_celltypes[num]
                        cmd += self.add_cmd_row(f'--celltype {white_celltype}')
                else:
                    exit("config.yaml 文件中 sub_clusters white_celltypes 填写格式错误 请查看注释信息！")
            if extraGene != 'None':
                cmd += self.add_cmd_row(f'-x {extraGene}')  
            if tissue != 'None':
                cmd += self.add_cmd_row(f'--tissue {tissue}')

            cmd += self.add_cmd_row(f'-o sub_{cell_name}/',end=True)
                
            with open(out_script,"w") as f:
                f.write(cmd)
            print(f"脚本 {out_script} 已生成") 

            num += 1
        #### 物种信息保存至数据库 
        db_update_bg = self.pjif  
        db_update_bg['species'] = species
        db_update_bg['tissue'] = tissue
        db_update_bg['custom_ref'] = custom_ref
        self.update_info_bag = db_update_bg
        
            
            
            
if __name__ == '__main__':
    pass
