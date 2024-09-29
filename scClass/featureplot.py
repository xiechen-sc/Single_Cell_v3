from .base_class import BaseClass
from single_cell_auto.util import *
import re
# featureplot
class Featureplot(BaseClass):
        analysis_module = 'featureplot'

        def get_script(self):
            outdir = self.outdir
            out_script = outdir + '/cmd_featureplot.sh'
            cpu = self.cpu
            plot = self.plot
            seurat = self.input_seurat
            assay = self.assay
            output = self.output
            genelist = self.genelist
            groupby = self.groupby
            reduct = self.reduct
            pvalue = self.pvalue
            splitby = self.splitby
            dotsplit = self.dotsplit
            selcet = self.selcet
            select_col = self.select_col
            select_lst = self.select_lst
            
            cmd = f"""set -e\nmodule purge && module load OESingleCell/3.0.d\n"""
            
            if selcet:
                select_lst_all = select_lst
            else:
                select_lst_all = ['all']
            for col_type in select_lst_all:
                if col_type != 'all':
                    output_dir = output + '_' + col_type
                else:
                    output_dir = output
                cmd += self.add_cmd_row('Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool')  
                cmd += self.add_cmd_row(f'-i {seurat}')
                cmd += self.add_cmd_row(f'-f h5seurat')
                cmd += self.add_cmd_row(f'-o {output_dir}')
                cmd += self.add_cmd_row(f'-j {cpu}') 
                cmd += self.add_cmd_row(f'--assay {assay}')
                cmd += self.add_cmd_row(f'--dataslot data')
                if type(col_type) == list:  # 多种细胞混合
                    sel_str = ""
                    for i in col_type:
                            sel_str += "\\'" + i + "\\',"
                    sel_str = re.sub(',$','',sel_str)                       
                    sel_col_str = "c(" + sel_str + ")"
                    cmd += self.add_cmd_row(f'--predicate "{select_col} %in% {sel_col_str}"')
                elif col_type == 'all':
                    pass
                else:
                    sel_col_str = "c(\\'" + col_type + "\\')"
                    cmd += self.add_cmd_row(f'--predicate "{select_col} %in% {sel_col_str}"')
                cmd += self.add_cmd_row(f'visualize')
                cmd += self.add_cmd_row(f'-x {genelist}') 
                cmd += self.add_cmd_row(f'-g {groupby}')
                cmd += self.add_cmd_row(f'--reduct {reduct}')
                cmd += self.add_cmd_row(f'-m {plot}')
                cmd += self.add_cmd_row(f'--vcolors customecol2')  
                cmd += self.add_cmd_row(f'--ccolors spectral')   
                cmd += self.add_cmd_row(f'--pointsize 0')
                if pvalue != 'None':
                    cmd += self.add_cmd_row(f'--pvalue {pvalue}') 
                if splitby != 'None':
                    cmd += self.add_cmd_row(f'--splitby {splitby}')
                cmd += self.add_cmd_row(f'--dotsplit {dotsplit}')  
                cmd += self.add_cmd_row(f'--dodge F',True)           

            f = open(out_script,'w')
            f.write(cmd)
            f.close()
            print(f"脚本 {out_script} 已生成")

