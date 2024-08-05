set -e
module purge && module load OESingleCell/3.0.d
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \
-i seurat.h5seurat \
-f h5seurat \
-o newcelltype \
-d h5seurat \
--update FALSE \
--assay RNA \
--dataslot counts,data,scale.data  \
changecelltype \
-c newcelltype.tsv \
-C clusters \
--palette customecol2 \
--reduct umap \
-b F
                
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \
-i  newcelltype/seurat.h5seurat \
-f h5seurat \
-o ./ \
--assay RNA \
--dataslot data \
summarize \
--reduct umap \
--palette customecol2 \
-c new_celltype \
-b sampleid,group \
--pointsize 0.5 \
--dosummary T

Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \
-i newcelltype/seurat.h5seurat \
-f h5seurat \
-o ./clusters_correlation \
-t 6 \
--assay RNA \
--slot data \
--reduct umap \
coefficient \
-g new_celltype

# 1.鉴定marker
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \
-i newcelltype/seurat.h5seurat \
-f h5seurat \
-o Marker \
--assay RNA \
--dataslot data,counts \
-j 10 \
findallmarkers \
-c 2 \
-N 10 \
-k 1 \
-p 0.05 \
-s F \
-e presto \
-n new_celltype

#2.可视化+anno
# marker热图
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis \
-i newcelltype/seurat.h5seurat \
-f h5seurat \
-o ./Marker \
-t 10 \
--assay RNA \
--slot data,scale.data \
heatmap \
-l Marker/top10_markers_for_each_cluster.xls \
-c gene_diff \
-n 10 \
-g new_celltype \
--group_colors customecol2 \
--sample_ratio 0.8 \
--style seurat

# featureplot小提琴图
Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool \
-i newcelltype/seurat.h5seurat  \
-f h5seurat \
-o ./Marker \
-j 10 --assay RNA \
--dataslot data \
visualize -l Marker/top10_markers_for_each_cluster.xls -g new_celltype \
--reduct umap \
--topn  10  \
--topby gene_diff \
-m vlnplot,featureplot \
--vcolors customecol2 \
--ccolors spectral \
--pointsize 0.3 \
--dodge F

#anno


Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \
-g Marker/all_markers_for_each_cluster.xls \
--anno /data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_annotation.xls  # 根据物种修改
Rscript  /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation \
-g Marker/top10_markers_for_each_cluster.xls \
--anno /data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_annotation.xls  # 根据物种修改
