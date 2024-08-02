from .base_class import BaseClass
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
            f = open(out_script,'w')
            f.write(f"""set -e
module purge && module load OESingleCell/3.0.d            
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \\
  -i {seurat}  \\
  -f h5seurat \\
  -o {output} \\
  -j {cpu} \\
  --assay {assay} \\
  --dataslot data \\
  visualize \\
  -x {genelist} \\
  -g {groupby} \\
  --reduct {reduct} \\
  -m {plot} \\
  --vcolors customecol2 \\
  --ccolors spectral \\
  --pointsize 0 \\
  --dodge F
""")
            f.close()
            print(f"脚本 {out_script} 已生成")

