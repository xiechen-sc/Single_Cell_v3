from .base_class import BaseClass
from single_cell_auto.util import cell_name_normalization,get_species_info,jinggao,database_add
from single_cell_auto.cmd_module import volcano
# diff
class Diff(BaseClass):

    def get_script(self):
            analysis_module = self.run
            input_rds = self.seurat
            cell_types = [str(i) for i in self.cell_types ]
            analysis_type = self.sub_type
            treat_all = self.treat
            control_all = self.control
            fc = self.fc
            p = self.p
            vs_type = self.vs_type
            species = self.species
            volcano_plot = self.volcano_plot

            top = self.top
            outdir = self.outdir

            len_t = len(treat_all); len_c = len(control_all)
            species_info = get_species_info(species=species)
            if len_t != len_c:
                exit('请检查 yaml 文件 实验组与对照组数量不一致')
            try:
                anno = species_info['anno']
            except KeyError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 的 anno 在数据库中不存在 请手动填写！')
            except TypeError:
                anno = '# 请手动填写！！！'
                jinggao(f'{species} 在数据库中不存在 请手动填写！')
            
            
            for num in range(0,len_t):
                cmd = "" # 定义一个空命令
                treat = treat_all[num]
                control = control_all[num]

                for cell_type in cell_types:
                    cell_type_out = cell_name_normalization(cell_type)
                    if cell_type != 'all':
                        if '[' in cell_type: # 列表
                            cell_type_out = cell_name_normalization(cell_type.replace("[","").replace("]","").replace(r",",r"_").replace(r" ",r"").replace(r"'",r""))
                        else:
                            cell_type_out = cell_name_normalization(cell_type)
                    cmd = f"""set -e
module purge
module load OESingleCell/3.0.d
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  \\
-i {input_rds}     \\
-f h5seurat     \\
-o ./{cell_type_out}-Diffexp/{treat}-vs-{control}     \\
--assay RNA     \\
--dataslot data,counts     \\
-j 10  \\
"""

                    if cell_type != 'all':
                        if '[' in cell_type: # 列表
                            sub_list = cell_type.replace("[","").replace("]","")
                            cell_type_out = cell_name_normalization(sub_list.replace(r",",r"_").replace(r" ",r"").replace(r"'",r""))
                            sub_list = sub_list.replace(r"'",r"\'")

                        else:
                            cell_type_out = cell_name_normalization(cell_type)
                            sub_list = "\\'" + cell_type + "\\'"

                        cmd = cmd + f"--predicate \"{analysis_type} %in% c({sub_list})\" \\\n"

                    cmd = cmd + f"""diffexp     \\
-c {vs_type}:{treat}:{control}     \\
-k {fc}     \\
-p {p}    \\
-e presto

Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  annotation \\
-g ./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-all_diffexp_genes.xls \\
--anno {anno}/annotation/gene_annotation.xls

Rscript   /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  annotation \\
-g ./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-diff-pval-{p}-FC-{fc}.xls \\
--anno {anno}/annotation/gene_annotation.xls

### diffexp_heatmap
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \\
-i {input_rds}  \\
-f h5seurat \\
-o ./{cell_type_out}-Diffexp/{treat}-vs-{control}  \\
-t 10 \\
--assay RNA \\
--slot data,scale.data \\
"""                 
                    if cell_type != 'all':
                        cmd += f"""--predicate "{analysis_type} %in% c({sub_list}) & {vs_type} %in% c(\\'{treat}\\',\\'{control}\\')" \\\n"""
                    cmd += f"""diff_heatmap \\
-d ./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-diff-pval-{p}-FC-{fc}.xls \\
-n {top} \\
-g {vs_type} \\
--group_colors customecol2 \\
--sample_ratio 0.8

rm ./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-all_diffexp_genes.xls ./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-diff-pval-{p}-FC-{fc}.xls



/public/scRNA_works/pipeline/scRNA-seq_further_analysis/enrichwrap.sh \\
-i {cell_type_out}-Diffexp/{treat}-vs-{control}/*-vs-*-diff-*.xls \\
-g  {anno} \\
-o {cell_type_out}-Diffexp/{treat}-vs-{control}/enrichment \\
-d TRUE
"""
                    
                    cmd_vol = volcano(
                        input=f'./{cell_type_out}-Diffexp/{treat}-vs-{control}/{vs_type}_{treat}-vs-{control}-all_diffexp_genes_anno.xls',
                        pvalue=p,
                        log2fc=fc,
                        output=f'./{cell_type_out}-Diffexp/{treat}-vs-{control}/'
                                )
                    if volcano_plot:
                        cmd += cmd_vol
                        
                    out_script = f'{outdir}/cmd_{cell_type_out}-{treat}-vs-{control}.diff.sh'
                    with open(out_script,"w") as f:
                        f.write(cmd)
                    print(f"脚本 {out_script} 已生成")
            #### 物种信息保存至数据库
            database_add(config_path=out_script,config_info={'species':species})
                    

