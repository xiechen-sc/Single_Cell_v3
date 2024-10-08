from .base_class import BaseClass
from single_cell_auto.util import *
from single_cell_auto.cmd_module import *
class Scenic(BaseClass):
    analysis_module = 'scenic'

    def get_step1_cmd(self,input_rds,db,scenic_species,coexMethod,output):
        cmd = '### scenic step 1 ###\n'
        cmd += self.add_cmd_row(f'Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/scenic.R')
        cmd += self.add_cmd_row(f'-i {input_rds}')
        cmd += self.add_cmd_row(f'-f seurat')
        cmd += self.add_cmd_row(f'-d {db}')
        cmd += self.add_cmd_row(f'-s {scenic_species}')
        cmd += self.add_cmd_row(f'--coexMethod {coexMethod}')
        cmd += self.add_cmd_row(f'-o {output}',end=True)
        return cmd
    def get_step2_cmd(self,step2_input,sub_seurat,seurat_sub_col,seurat_sub_col_value,result_dir,rss_rank_top_gene,groupby):
        cmd = '\n### scenic step 2 ###\n'
        lst = ['all']  # 随机值 不用管
        if sub_seurat: 
            lst = seurat_sub_col_value
        # ras rss csi
        for sub_value in lst:
            s_v = sub_value
            cmd_predicate = ''
            if s_v == 'all':
                out_dir = f'scenic_ras_rss_by_{groupby}'
                out_dir_csi = f'scenic_csi_by_{groupby}'
            elif type(s_v) == list:
                predicate = "\'" + "\',\'".join(s_v) + "\'"
                outname = '_'.join(s_v)
                cmd_predicate = f'--predicate "{seurat_sub_col} %in% c({predicate})'
                out_dir = f'scenic_ras_rss_{seurat_sub_col}_{outname}_by_{groupby}'
                out_dir_csi = f'scenic_csi_{seurat_sub_col}_{outname}_by_{groupby}'
            else:
                predicate = '\'' + s_v + '\''
                cmd_predicate = f'--predicate "{seurat_sub_col} %in% c({predicate})'
                out_dir = f'scenic_ras_rss_{seurat_sub_col}_{s_v}_by_{groupby}'
                out_dir_csi = f'scenic_csi_{seurat_sub_col}_{s_v}_by_{groupby}'
            cmd_ras_rss = ''
            cmd_ras_rss += self.add_cmd_row(f'Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/RunRAS-RSS.R')
            cmd_ras_rss += self.add_cmd_row(f'-i {step2_input}')
            cmd_ras_rss += self.add_cmd_row(f'-v {result_dir}')
            cmd_ras_rss += self.add_cmd_row(f'-f rds')
            cmd_ras_rss += self.add_cmd_row(f'-t {rss_rank_top_gene}')
            cmd_ras_rss += self.add_cmd_row(f'-c {groupby}')
            cmd_ras_rss += self.add_cmd_row(f'-s 0')
            cmd_csi = ''
            cmd_csi += self.add_cmd_row(f'Rscript /home/luyao/10X_scRNAseq_v3/src/GRN/RunCSI.R')
            cmd_csi += self.add_cmd_row(f'-i {step2_input}')
            cmd_csi += self.add_cmd_row(f'-v {result_dir}')
            cmd_csi += self.add_cmd_row(f'-f rds')
            cmd_csi += self.add_cmd_row(f'-c {groupby}')
            cmd_csi += self.add_cmd_row(f'-n {self.cluster_n}')
            if cmd_predicate:
                cmd_ras_rss += self.add_cmd_row(cmd_predicate)
                cmd_csi += self.add_cmd_row(cmd_predicate)
            if self.use_color_anno:
                cmd_ras_rss += self.add_cmd_row(f'--use_color_anno {self.use_color_anno}')
            if self.color_file:
                cmd_ras_rss += self.add_cmd_row(f'--color_file {self.color_file}')
            if self.palette:
                cmd_ras_rss += self.add_cmd_row(f'--palette  {self.palette}')
            cmd_ras_rss += self.add_cmd_row(f'-o {out_dir}',end=True)
            cmd_csi += self.add_cmd_row(f'-o {out_dir_csi}',end=True)
            cmd += cmd_ras_rss
            cmd += cmd_csi

        return cmd
            
                


        


    def get_script(self):
        step1_run = self.step1_run
        step1_input = self.step1_input
        species = self.species
        coexMethod = self.coexMethod
        step1_outdir = self.step1_outdir
        step2_run = self.step2_run
        step2_input = self.step2_input
        result_dir = self.result_dir
        rss_rank_top_gene = self.rss_rank_top_gene
        groupby = self.groupby
        sub_seurat = self.sub_seurat
        seurat_sub_col = self.seurat_sub_col
        seurat_sub_col_value = self.seurat_sub_col_value
        use_color_anno = self.use_color_anno
        color_file = self.color_file
        palette = self.palette
        cluster_n = self.cluster_n
        analysis_type = self.analysis_module
        # species 处理
        if species not in Frequent_species():
            jinggao('物种只能是人或小鼠！')
            exit()
        else:
            if 'human' in species:
                scenic_species = 'hgnc'
            elif 'mouse' in species:
                scenic_species = 'mgi'
            db_file = get_species_info(species=species)
            db = db_file['scenic']
        cmd = 'set -e\nmodule purge\nsource /home/lipeng/miniconda3/bin/activate Scenic\n'
        cmd_seurat2rds,step1_input_2 = self.judgment_seurat_rds(input_file=step1_input)  # step1 转化为了 rds
        if step1_run == False and step2_run == False:
            jinggao("step1 step2 至少执行一个！！！")
            exit()
        if cmd_seurat2rds != '':
            cmd += cmd_seurat2rds
        if step1_run:
            cmd += self.get_step1_cmd(input_rds=step1_input_2,db=db,scenic_species=scenic_species,coexMethod=coexMethod,output=step1_outdir)  # 这个 step1_outdir 需要添加至数据库 添加绝对路径
        if step2_run:
            if step1_input == step2_input:
                step2_input_2 = step1_input_2
            else:
                cmd_seurat2rds_step2,step2_input_2 = self.judgment_seurat_rds(input_file=step2_input)  # step1 转化为了 rds
                cmd += cmd_seurat2rds_step2
            cmd += self.get_step2_cmd(step2_input=step2_input_2,sub_seurat=sub_seurat,seurat_sub_col=seurat_sub_col,seurat_sub_col_value=seurat_sub_col_value,result_dir=result_dir,rss_rank_top_gene=rss_rank_top_gene,groupby=groupby)
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        with open(out_script,"w") as f:
                f.write(cmd)
        print(f"脚本 {out_script} 已生成")
        # 更新数据库  
        import os
        db_update_bg = self.pjif  # 可能有 scenic 也可能没有
        if 'scenic' not in db_update_bg:
            db_update_bg['scenic'] = dict()
        if step1_run:
            rst_path = os.path.abspath(step1_outdir)
            rst_path += '/int/3.4_regulonAUC.Rds'
            if 'step1_rst' not in db_update_bg['scenic']:
                db_update_bg['scenic']['step1_rst'] = [rst_path]
            else:
                db_update_bg['scenic']['step1_rst'].append(rst_path)
            db_update_bg['scenic']['scenic_species'] = scenic_species
            db_update_bg['scenic']['step1_rds'] = os.path.abspath(step1_input_2) 
        
        self.update_info_bag = db_update_bg

        

