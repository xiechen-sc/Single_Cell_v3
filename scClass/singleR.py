from .base_class import BaseClass
from single_cell_auto.util import *
from single_cell_auto.cmd_module import *
import re
class SingleR(BaseClass):

    def get_script(self):
        outdir = self.outdir
        seurat = self.seurat
        output = self.output
        prefix = self.result_perfix
        assay = self.assay
        singleR_rds = self.singleR_rds
        reduct2 = self.reduct2
        species = self.species
        species_input2singleR = re.sub('_.*','',species)  
        annolevel = self.annolevel
        species_info = get_species_info(species=species)
        if singleR_rds == 'default':
            try:
                singleR_rds = species_info['singleR'][singleR_rds]
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
        prefix_num = len(prefix)
        seurat_num = len(seurat)
        
        if prefix_num == 0 and seurat_num == 1:
            out_script = outdir + '/cmd_singleR.sh'
            seurat = seurat[0]
            cmd = cmd_singleR(seurat=seurat,assay=assay,output=output,singleR_rds=singleR_rds,reduct2=reduct2,species=species_input2singleR,annolevel=annolevel)
            f = open(out_script,'w')
            f.write(cmd)
            f.close()
            print(f"脚本 {out_script} 已生成")
        elif prefix_num != seurat_num:
            jinggao('config 中 输入多个 seurat 时，必须有相同数量的 prefix 对应')
            exit(1)
        else:  # 多个 seurat 与 多个 prefix
            for i in range(0,seurat_num):
                seurat_1 = seurat[i]
                prefix_1 = prefix[i]
                output_1 = output + f'/{prefix_1}'
                cmd = cmd_singleR(seurat=seurat_1,assay=assay,output=output_1,singleR_rds=singleR_rds,reduct2=reduct2,species=species_input2singleR,annolevel=annolevel)
                out_script = outdir + f'/cmd_singleR_{prefix_1}.sh'
                f = open(out_script,'w')
                f.write(cmd)
                f.close()
                print(f"脚本 {out_script} 已生成")
                

