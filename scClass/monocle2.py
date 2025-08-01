from .base_class import BaseClass
from single_cell_auto.util import *

class Monocle2(BaseClass):
    analysis_module = "monocle2"

    def get_step1_cmd(self,input_rds,assay,col_name,sub_seurat,sub_col,sub_lst,step1_groupby,output,cores_use,resolution,result_rds):
        cmd = '### run monocle2 ###\n'
        cmd += self.add_cmd_row(f'Rscript /home/luyao/10X_scRNAseq_v3/src/Pseuduotime/monocle.R')
        cmd += self.add_cmd_row(f'-i {input_rds}')
        cmd += self.add_cmd_row(f'-f seurat')
        cmd += self.add_cmd_row(f'--assay {assay}')
        cmd += self.add_cmd_row(f'-d {col_name}')
        cmd += self.add_cmd_row(f'-C {step1_groupby}')
        if sub_seurat: 
            str_list = [str(k) for k in sub_lst.split(",")]
            cell_name = "_".join(str_list)
            cell_type = ",".join(["'" + k + "'" for k in str_list])
            cell_name_raw = str_list
            cmd += self.add_cmd_row(f'--predicate  "{sub_col} %in% c({cell_type})"')
            # cmd += self.add_cmd_row(f'-c {sub_col}')
            # cmd += self.add_cmd_row(f'-u {sub_lst}')
        cmd += self.add_cmd_row(f'-j {cores_use} -x 0.01 -r {resolution} -s 1')
        if result_rds:
            cmd += self.add_cmd_row(f'--rds {result_rds}')
        if self.use_color_anno:
                cmd += self.add_cmd_row(f'--use_color_anno TRUE')
        if self.color_file:
                cmd += self.add_cmd_row(f'--color_file {self.color_file}')
        if self.palette:
                cmd += self.add_cmd_row(f'--palette  {self.palette}')
        cmd += self.add_cmd_row(f'-o {output}',end=True)
        return cmd

    def get_step2_cmd(self,monocle_rds,genelist,root_state,show_branch,branch,vis_methods,step2_groupby,anno,output,module_expressplot):
        cmd = '\n### monocle2 下游 ###\n'
        cmd += self.add_cmd_row(f'Rscript  /public/scRNA_works/pipeline/scRNA-seq_further_analysis/visualize_pseudotime.R')
        cmd += self.add_cmd_row(f'-i {monocle_rds}')
        cmd += self.add_cmd_row(f'-g {genelist}')
        cmd += self.add_cmd_row(f'--root {root_state}')
        if show_branch:
            cmd += self.add_cmd_row(f'--show_branch {show_branch}')
        if branch:
            cmd += self.add_cmd_row(f'-b {branch}')
        cmd += self.add_cmd_row(f'-m {vis_methods}')
        cmd += self.add_cmd_row(f'-c {step2_groupby} -j 5')
        cmd += self.add_cmd_row(f'-a {anno}/annotation/gene_annotation.xls')
        if module_expressplot:
            cmd += self.add_cmd_row(f'--module {module_expressplot}')
        if self.use_color_anno:
                cmd += self.add_cmd_row(f'--use_color_anno TRUE')
        if self.color_file:
                cmd += self.add_cmd_row(f'--color_file {self.color_file}')
        if self.palette:
                cmd += self.add_cmd_row(f'--palette  {self.palette}')
        cmd += self.add_cmd_row(f'-o {output}',end=True)
        return cmd

    def get_enrichment_cmd(self,output,anno):
        cmd = '\n### 各module基因富集分析 ###\nmodule purge\nmodule load OESingleCell/3.0.d\n'
        cmd += self.add_cmd_row(f'/public/scRNA_works/pipeline/scRNA-seq_further_analysis/enrichwrap.sh')
        cmd += self.add_cmd_row(f'-i {output}/pseudotime_heatmap_gene_module_anno.xls')
        cmd += self.add_cmd_row(f'-g {anno}')
        cmd += self.add_cmd_row(f'-o {output}/enrich')
        cmd += self.add_cmd_row(f'-p TRUE -n 2',end=True)
        return cmd

    def get_script(self):
        monocle_run = self.step1_monocle_run
        input_seurat = self.input_seurat
        assay = self.assay
        col_name = self.col_name
        step1_groupby = self.step1_groupby
        sub_seurat = self.sub_seurat
        sub_col = self.sub_col
        sub_lst = self.sub_lst
        cores_use = self.cores_use
        resolution = self.resolution
        output_dir = self.output_dir
        result_rds = self.result_rds
        step2_downstream_run = self.step2_downstream_run
        monocle_rds = self.monocle_rds
        genelist = self.genelist
        root_state = self.root_state
        show_branch =self.show_branch
        branch = self.branch
        vis_methods = self.vis_methods
        step2_groupby = self.step2_groupby
        module_enrichment = self.module_enrichment
        module_expressplot = self.module_expressplot
        species = self.species
        species_info = get_species_info(species=species)
        try:
            anno = species_info['anno']
        except KeyError:
            anno = '# 请手动填写！！！'
            jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
        except TypeError:
            anno = '# 请手动填写！！！'
            jinggao(f'{species} 在数据库中不存在 请手动填写！')
        cmd = 'set -e\nmodule purge\nmodule load OESingleCell/2.0.0\n'
        if monocle_run:
            cmd_seurat2rds,input_rds = self.judgment_seurat_rds(input_file=input_seurat)  # 如果输入文件为h5seurat,转化为rds
            if cmd_seurat2rds != '':
                cmd += cmd_seurat2rds
            cmd += self.get_step1_cmd(input_rds=input_rds,assay=assay,col_name=col_name,sub_seurat=sub_seurat,sub_col=sub_col,sub_lst=sub_lst,step1_groupby=step1_groupby,output=output_dir,cores_use=cores_use,resolution=resolution,result_rds=result_rds)
        if step2_downstream_run:
            cmd += self.get_step2_cmd(monocle_rds=monocle_rds,genelist=genelist,root_state=root_state,show_branch=show_branch,branch=branch,vis_methods=vis_methods,step2_groupby=step2_groupby,anno=anno,output=output_dir,module_expressplot=module_expressplot)
        if module_enrichment:
            cmd += self.get_enrichment_cmd(output=output_dir,anno=anno)
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        with open(out_script,"w") as f:
                f.write(cmd)
        print(f"脚本 {out_script} 已生成")



        

