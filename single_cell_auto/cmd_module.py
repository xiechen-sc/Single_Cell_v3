from single_cell_auto.util import jinggao,get_species_info
# singleR 注释
def cmd_singleR(seurat,output,assay,singleR_rds,reduct2,species,annolevel):
    cmd = f"""set -e
module purge && module load OESingleCell/3.0.d
# reference celltype
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool  \\
-i {seurat}  \\
-f h5seurat \\
-o {output} \\
-d h5seurat \\
--update T \\
--assay {assay} \\
--dataslot counts \\
celltyping \\
-r {singleR_rds} \\
--annolevel {annolevel} \\
--usecluster F \\
--demethod classic \\
--pointsize 0.3 \\
-n 25 \\
--reduct {reduct2} \\
--species {species}
"""
    return cmd


def volcano(input,sig,pvalue,log2fc,output,symbol_topn):
    if sig == 'qval':
        psub = f"-q {pvalue}"
    elif sig == 'pval':
        psub = f"-p {pvalue}"
    if symbol_topn:
        cmd = f"""
set -e
module purge && module load OESingleCell/3.0.d
Rscript /gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/volcanoplot/volcano.r \\
-i {input} \\
{psub} \\
-f {log2fc} \\
-o {output} 
"""
    else:
        cmd = f"""
set -e
module purge && module load OESingleCell/3.0.d
Rscript /gpfs/oe-scrna/pipeline/scRNA-seq_further_analysis/volcanoplot/volcano.r \\
-i {input} \\
{psub} \\
-f {log2fc} \\
-o {output} \\
--symbol_topn {symbol_topn}
        """
    return cmd


    
    