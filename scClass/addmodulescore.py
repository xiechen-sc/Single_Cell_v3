from .base_class import BaseClass
from single_cell_auto.util import cell_name_normalization,get_species_info,jinggao,database_add
# Addmodulescore
class Addmodulescore(BaseClass):

    def get_script(self):
            analysis_module = self.run
            outdir = self.outdir
            input = self.input # 输入 seurat 对象文件
            species = self.species  # 物种
            genelist = self.genelist  # 输入genelist 列表
            groupby = self.groupby  # 分组方式
            splitby = self.splitby  # 拆分展示小提琴图
            fsplitby = self.fsplitby  # 拆分展示umap图
            pvalue = self.pvalue
            sub = self.sub
            
            # 下方内容选择性填写 建议默认 
            show_box = self.show_box
            reduct = self.reduct # 降维方式 
            output = self.output  # 输出目录
            assay = self.reduct  # RNA  SCT
            dataslot = self.dataslot  # counts,data,scale.data
            strict = self.strict  # 是否使用严格模式筛选gene，默认FALSE
            assay = self.assay
            pointsize = self.pointsize
            species_info = get_species_info(species=species)

            try:
                anno = species_info['anno']
            except KeyError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
            except TypeError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 在数据库中不存在 请手动填写！')

            if sub == 'None':
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
                    cmd += f'    --splitby {splitby} \\\n'
                if not (fsplitby == 'None'):
                    cmd += f'    --fsplitby {fsplitby} \\\n'
                if show_box:
                    cmd += f'    --show_box {show_box} \\\n'
                if not (strict == 'None'):
                    cmd += f'    	--strict {strict} \\\n'
                if not (pvalue == 'None'):
                    cmd += f'    	--pvalue {pvalue} \\\n'
                cmd += f'        --pointsize {pointsize} '

                out_script = f'{outdir}/cmd_addmodulescore.sh'
                with open(out_script,"w") as f:
                    f.write(cmd)
                print(f"脚本 {out_script} 已生成")
            else:
                select =  self.select
                groups_sel = []
                for sel in select:
                    if isinstance(sel,list):
                        group = '\',\''.join(map(str, sel))
                        group = '\'' + group +'\''
                    else:
                        group = '\'' + str(sel) +'\''
                    
                    groups_sel.append(group)
                for group in groups_sel:
                    import re
                    group2 = re.sub(r'\'|\"','',group)
                    group2 = re.sub(r'\,','_',group2)
                    output2 = f"{sub}_{group2}_addmodulescore"
                    group3 = re.sub(r'\'','\\\'',group)
                    cmd = f"set -e\nmodule purge && module load OESingleCell/3.0.d\n"
                    cmd += f"""Rscript /gpfs/oe-scrna/chenhaoruo/script/new-sctool/sctool \\
    --input {input} \\
    --reduct {reduct} \\
    --output {output2} \\
    --assay {assay} \\
    --dataslot {dataslot} \\
    --anno {anno} \\
    --predicate \"{sub} %in% c({group3})" \\
    sc_addmodulescore  \\
        -x {genelist} \\
        -g {groupby} \\
"""
                    if not (splitby == 'None'):
                        cmd += f'    	--splitby {splitby} \\\n'
                    if not (fsplitby == 'None'):
                        cmd += f'    	--fsplitby {fsplitby} \\\n'
                    if show_box:
                        cmd += f'    	--show_box {show_box} \\\n'
                    if not (strict == 'None'):
                        cmd += f'    	--strict {strict} \\\n'
                    if not (pvalue == 'None'):
                        cmd += f'    	--pvalue {pvalue} \\\n'
                    cmd += f'        --pointsize {pointsize} '
                    
                    out_script = f'{outdir}/cmd_{sub}_{group2}_addmodulescore.sh'
                    with open(out_script,"w") as f:
                        f.write(cmd)
                        print(f"脚本 {out_script} 已生成")

                    


            #### 物种信息保存至数据库
            db_update_bg = self.pjif
            db_update_bg['species'] = self.species
            self.update_info_bag = db_update_bg
                        

