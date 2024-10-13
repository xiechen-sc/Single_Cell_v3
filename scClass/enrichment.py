from .base_class import BaseClass
from single_cell_auto.util import *
import os
class Enrichment(BaseClass):
        analysis_module = 'enrichment'

        def get_top_marker_file(self,top_n,input,outdir,prefix,sort_by):
            import pandas as pd
            df = pd.read_csv(input,sep='\t')
            cell_type = df.columns.to_list()[6]
            print(f'正在对 {input} 执行提取 top {top_n} 处理 请等待！')
            df_sorted = df.reindex(df[sort_by].abs().sort_values().index)  # 根据 绝对值大小进行排序
            df_top_n = df_sorted.groupby(cell_type).head(top_n)
            ipt = outdir + f'/{prefix}_all_clusters_marker_top_{top_n}.xls'
            df_top_n.to_csv(ipt, index=False,sep='\t')
            print(f'{input} 提取 top {top_n} 后的文件为 {ipt} 脚本已更新脚本 input 文件！')
            return ipt
        def get_cmd(self,analysis_model,input,anno,prefix,top_n):
            cmd = f"""\n"""
            cmd += self.add_cmd_row(f"/public/scRNA_works/pipeline/scRNA-seq_further_analysis/enrichwrap.sh")
            cmd += self.add_cmd_row(f"-i {input}")
            cmd += self.add_cmd_row(f"-g {anno}")
            if analysis_model == '0':
                out_dir = prefix + '_genelist' 
                cmd += self.add_cmd_row(f"-e TRUE")
            elif analysis_model == '1':
                out_dir = prefix + f'_top_{top_n}' 
                cmd += self.add_cmd_row(f"-m TRUE")
                cmd += self.add_cmd_row(f"-n 7")
            elif analysis_model == '2':
                out_dir = prefix + '_diffexp_list' 
                cmd += self.add_cmd_row(f"-d TRUE")  
            elif analysis_model == '3':
                out_dir = prefix + '_diffexp_dir' 
            elif analysis_model == '4':
                out_dir = prefix + '_monocle2' 
                cmd += self.add_cmd_row(f"-p TRUE")  
                cmd += self.add_cmd_row(f"-n 2")  
            else:
                jinggao('出现了未知的分析类型 请检查 config 中的 analysis_model 参数 只能使用 0 1 2 3 4')
                exit(1)
            cmd += self.add_cmd_row(f"-o {out_dir}",end=True)  
            return cmd

        def get_script(self):
            input = self.input
            prefix_lst = self.prefix_lst
            sort_by = self.sort_by
            analysis_model = self.analysis_model
            species = self.species
            outdir = self.outdir
            top_n = self.top_n
            out_script = outdir + '/cmd_enrichment.sh'
            len_analysis_model = len(analysis_model)
            len_input = len(input)
            len_prefix_lst = len(prefix_lst)
            if len_prefix_lst != len_input and len_prefix_lst != 0:
                jinggao('config.yaml 中 prefix_lst 列表中数量必须与 input列表中元素个数一致或者 或者不填！')
            if len_analysis_model == 1:
                analysis_model = [str(analysis_model[0]) for i in range(len(input))]
            else:
                if len_input != len_analysis_model:
                    jinggao('config.yaml 中 analysis_model 列表中数量必须与 input列表中元素个数一致或者 数量为1！')
                    exit(1)
                analysis_model = [str(i) for i in analysis_model]
            species_info = get_species_info(species=species)
            try:
                anno = species_info['anno']
            except KeyError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
            except TypeError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 在数据库中不存在 请手动填写！')
            cmd = f"""set -e\nmodule purge && module load OESingleCell/3.0.d\n"""
            for i in range(len_input):
                input_i = input[i]
                analysis_model_i = analysis_model[i]
                if len_prefix_lst == 0:
                    prefix = 'all'
                else:
                    prefix = str(prefix_lst[i])


                if analysis_model_i == '1':
                    input_i = self.get_top_marker_file(input=input_i,outdir=outdir,prefix=prefix,top_n=top_n,sort_by=sort_by)
                cmd += self.get_cmd(analysis_model=analysis_model_i,input=input_i,anno=anno,prefix=prefix,top_n=top_n)
            with open(out_script,"w") as f:
                f.write(cmd)
            #### 物种信息保存至数据库
            db_update_bg = self.pjif  
            db_update_bg['species'] = species
            self.update_info_bag = db_update_bg
            print(f"脚本 {out_script} 已生成")