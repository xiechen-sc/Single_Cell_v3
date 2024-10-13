from .base_class import BaseClass
from single_cell_auto.util import cell_name_normalization,get_species_info,jinggao,database_add

# diff
class DecontX(BaseClass):
    analysis_module = 'decontX'
    def get_script(self):
        input_seurat = self.input_seurat
        threshold = self.threshold
        out_script = f'{self.outdir}/cmd_{self.analysis_module}.sh'
        reduct = self.reduct

        cmd = 'set -e\n'
        cmd += self.add_cmd_row('/gpfs/oe-scrna/zhengfuxing/pipline_self/decontX/decontX.r')
        
        if threshold != None:
            cmd += self.add_cmd_row(f'-t {threshold}')
        if reduct != None:
            cmd += self.add_cmd_row(f'-r {threshold}')
        cmd += self.add_cmd_row(f'-i {input_seurat}',end=True)
        with open(out_script,"w") as f:
                f.write(cmd)
        print(f"脚本 {out_script} 已生成")
        