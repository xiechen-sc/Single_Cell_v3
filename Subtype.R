#!/usr/bin/env Rscript 

# the command line api of data analysis of biological experiments
if (!requireNamespace('argparse', quietly = TRUE)) {
  stop("Please install argparse to parse the command line paramters!")
}
suppressPackageStartupMessages( library("argparse") )
suppressPackageStartupMessages( library("magrittr") )
suppressPackageStartupMessages( library("stringr") )

# ****************** COMMAND LINE PARAMETERS SETTING ****************
# ======================= GLOBAL command line parameters setting ======================= #
parser = ArgumentParser(description = "single cell analysis results visualization toolsets.", usage = "%(prog)s [global options]" )
parser$add_argument("-i", "--input", 
                    type = "character",
             help = "输入的需要进行分析的rds或者h5seurat文件")
parser$add_argument("-f", "--informat", 
                    type = "character", default = "h5seurat",
             help = "输入文件格式:h5seurat,或者rds.[default: %(default)s]")
parser$add_argument("-o", "--output",
                    type = "character", default = "./",
             help = "the output directory of results."  )
parser$add_argument("--assay", 
                    type = "character", default = "RNA",
             help = "分析用的assay.[default %(default)s]")
parser$add_argument("-m", "--method", 
                    type = "character", default = "clustering,marker,celltyping",
             help = "做亚群分析的步骤，逗号分隔，可选clustering,marker,celltyping,marker_vis.[default %(default)s]")
parser$add_argument("--predicate", 
                    type = "character", default = NULL,
             help = "The conditional expression to subset cells used for subtyping.[default: %(default)s]")
parser$add_argument("--subset", 
                    type = "character", default = NULL,
             help = "The conditional expression to subset cells used for subtyping.[default: %(default)s]")
parser$add_argument("-s", "--species", 
                    type="character", default = NULL,
             help="分析的物种，可分为human,mouse,rat,other.[default: %(default)s]")
parser$add_argument("--pointsize", 
                    type="character", default = "0.5",
             help="聚类图散点大小.[default: %(default)s]")
parser$add_argument("--reduct1", 
                    type = "character", default = "pca",
             help = "[OPTIONAL]亚群分析降维聚类的一级降维方式.[default %(default)s]")
parser$add_argument("--reduct2",
                    type = "character", default = "umap",
             help = "[OPTIONAL]亚群分析降维聚类的二级降维方式.[default %(default)s]")
parser$add_argument("--resolution",
                    type = "character", default = "0.4",
             help = "[OPTIONAL]亚群分析降维聚类的分辨率.[default %(default)s]")
parser$add_argument("--batchid", 
                    type = "character", default = "batchid",
             help = "[OPTIONAL]亚群分析降维聚类的一级降维去批次时候用的批次信息.[default %(default)s]")
parser$add_argument("--anno", 
                    type = "character", default = NULL,
             help = "[OPTIONAL]亚群分析物种的基因组信息.[default %(default)s]")
parser$add_argument("-r","--customref", 
                    type = "character",default = NULL,
             help = "[OPTIONAL]亚群分析细胞类型鉴定所需的参考文件，若无则根据物种进行默认选择.[default %(default)s].")
parser$add_argument("--annolevel", 
                    type = "character",default = "single",
             help = "[OPTIONAL]亚群分析细胞类型鉴定参考文件的细胞类型精细度.[default %(default)s].")
parser$add_argument("--celltype", 
                    type="character", default = NULL,
             help = "手动指定用哪个细胞类型的经典marker出基因可视化结果，逗号分隔，需要一字不差完全匹配，T_cell,B_cell,NK,Endothelial,Epithelial,Fibroblast,Macrophage,Microglia,Neuron,Neutrophil,Myeloid.[default: %(default)s]")
parser$add_argument("--tissue", 
                    type="character", default = NULL,
             help = "指定物种的组织器官类型，可选择(英文填写)：brain(脑)、Intestinal(肠)、lung(肺)、gastric(胃癌)、tumour(肿瘤).[default: %(default)s]")
parser$add_argument("--default_markers", 
                    type="character", default = "TRUE",
             help = "是否要出默认marker基因列表的基因可视化.[default: %(default)s]")
parser$add_argument("-x", "--extraGene", 
                    type = "character", default = NULL,
             help = "[OPTIONAL]除了经典marker基因，还可以输入额外的基因文件，文件内容tab分隔，多个marker基因文件可用逗号分隔传参.")

opt = parser$parse_args()
# ======================= GLOBAL command line parameters setting ======================= #
# 做亚型筛选
if ( !is.null(opt$subset) ){
    condition = unlist(strsplit(opt$subset, ";"))
    all_types = unlist(lapply(condition, function(x){ 
                             a = unlist(strsplit(x, ":"))[2]
                             b = unlist(strsplit(a, ","))
                             b = gsub("[^[:alnum:]_]", "_", b)
                             return(b)}))
    filename = paste(all_types, collapse="_")
    all_condition = unlist(lapply(condition, function(x){ 
                             a = unlist(strsplit(x, ":"))
                             b = unlist(strsplit(a[2], ","))
                             elements_str = paste0("\\\'", b, "\\\'")
                             result = paste0(a[1]," %in% c( ",paste(elements_str, collapse = ", "),")" )
                             return(result)}))
    predicate = paste0(all_condition, collapse = " & ")
    predicate = paste0("--predicate ", shQuote(predicate, type = "sh") )
} else if ( !is.null(opt$predicate) ){
    pattern <- "(?<=%in% |== ).+?(?=\\s|\\||&|$)"
    types <- unlist(str_extract_all(opt$predicate, pattern))
    all_types = unlist(lapply(types, function(x){ a = eval(parse(text = x)); return(a)}))
    filename = paste(all_types, collapse="_")
    predicate = gsub("[\"|']","\\\\'", opt$predicate)
    predicate = paste0("--predicate ", shQuote(predicate, type = "sh") )
} else { 
    filename = "seurat"
    predicate = ""
}
# 亚型分析的方法
if ( is.null(opt$method) ){
    print("未提供有效的亚群分析方法，默认进行clustering、marker、celltyping分析")
    method = c("clustering","marker","celltyping")
}else if( opt$method == "all" ){
    method = c("clustering","marker","celltyping","marker_vis")
}else{
    method = unlist(strsplit(opt$method,","))
}

# output directory setting
if ( is.null(opt$output) ){
    print("NO output directory specified,the current directory will be used!")
    root_dir = getwd()
}else{
    if ( file.exists(opt$output)){
        root_dir = opt$output
    }else{
        root_dir = opt$output
        dir.create(root_dir, recursive = T)
    }
}
root_dir = normalizePath(root_dir)
#=================================================================================
# 组合运行命令
#=================================================================================
cmd = NULL
## 降维聚类、分样本可视化、细胞群相关性分析的命令组合
if ( "clustering" %in% method ){
    sub_h5seurat = paste0(root_dir,"/Clustering/",filename,".h5seurat")
    if ( opt$reduct1 == "mnn") { 
        reduct1 = paste0( opt$reduct1, " --batchid ",opt$batchid, " --components 10 ")
    } else if ( opt$reduct1 == "pca,harmony" ){ 
        reduct1 = paste0( opt$reduct1, " --batchid ",opt$batchid, " -t 20 -y 30 ")
    } else {
        reduct1 = opt$reduct1 
    }
    param = paste0("-i ",opt$input," -f ",opt$informat,
                   " -o ",root_dir,"/Clustering  -d h5seurat --assay ",opt$assay, " --prefix ",filename,
                   "  --dataslot counts,data,scale.data ",predicate," --update F  bclust --reduct1 ",reduct1," --reduct2 ",opt$reduct2,
                   " --clusteringuse snn  --resolution ",opt$resolution," --rerun T   --pointsize ",opt$pointsize," --palette customecol2")

    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", param) 
    cmd = c(cmd,"# 降维聚类", cmd1)

    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", 
                        " -i ",sub_h5seurat," -f h5seurat ",
                        " -o ",root_dir,"/Clustering  --assay ",opt$assay, " --dataslot data", 
                        " summarize --reduct ",opt$reduct2," --palette customecol2 -c clusters  -b sampleid,group --pointsize ",opt$pointsize," --dosummary T")
    cmd = c(cmd, "# 分样本可视化", cmd1)
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis ", 
                        " -i ",sub_h5seurat," -f h5seurat -t 6 ",
                        " -o ",root_dir,"/Clustering/clusters_correlation  --assay ",opt$assay, " --slot data", 
                        " coefficient  -g clusters")
    cmd = c(cmd, "# 细胞群相关性分析", cmd1)
    predicate = ""
} else { sub_h5seurat = opt$input }
## marker分析的命令组合
if ( "marker" %in% method ){
    param = paste0("-i ",sub_h5seurat," -f h5seurat",
                   " -o ",root_dir,"/Marker  -d h5seurat   --assay ",opt$assay, 
                   "  --dataslot data,counts ",predicate," findallmarkers -c 2  -N 10 -k 1 -p 0.05 -s F -e presto -n clusters")
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", param)
    cmd = c(cmd, "# marker基因鉴定", cmd1)

    param = paste0(" -i ",sub_h5seurat," -f h5seurat",
                   " -o ",root_dir,"/Marker   -t 10  --assay ",opt$assay, 
                   "  --slot data,scale.data ",predicate," heatmap -l ",root_dir,"/Marker/top10_markers_for_each_cluster.xls ",
                   " -c gene_diff -n 10 -g clusters --group_colors customecol2 --sample_ratio 0.8 --style seurat" )
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/scVis ", param)
    cmd = c(cmd, "# top marker基因热图", cmd1)

    param = paste0("-i ",sub_h5seurat," -f h5seurat",
                   " -o ",root_dir,"/Marker -j 10   --assay ",opt$assay, 
                   "  --dataslot data ",predicate," visualize -l ",root_dir,"/Marker/top10_markers_for_each_cluster.xls ",
                   " -g clusters --reduct ",opt$reduct2," --topn 10 --topby gene_diff -m vlnplot,featureplot",
                   " --vcolors customecol2 --ccolors spectral  --pointsize ",opt$pointsize," --dodge F")
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", param)
    cmd = c(cmd, "# top marker基因featureplot和vlnplot", cmd1)
    # anno
    if ( !tolower(opt$species) %in% c("human", "mouse") & is.null(opt$anno) ){
        print("特殊物种，未提供基因组注释文件，不进行基因注释!")
    } else {
        if ( !is.null(opt$anno)){
            anno = opt$anno
        } else if ( tolower(opt$species) %in% c("human", "mouse") & is.null(opt$anno) ){
            print("人或小鼠未指定基因组注释文件，默认使用2020版本的基因组信息!")
            anno = ifelse(tolower(opt$species) == "human","/data/database/cellranger-refdata/refdata-gex-GRCh38-2020-A/annotation/gene_annotation.xls",
                                                          "/data/database/cellranger-refdata/refdata-gex-mm10-2020-A/annotation/gene_annotation.xls")
        }
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation -g ",root_dir,"/Marker/all_markers_for_each_cluster.xls --anno ",anno)
    cmd = c(cmd, "# marker基因注释", cmd1)
    cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool annotation -g ",root_dir,"/Marker/top10_markers_for_each_cluster.xls --anno ",anno)
    cmd = c(cmd, cmd1)
    }
}
## SingleR细胞类型鉴定的命令组合
if ( "celltyping" %in% method ){
    if ( !tolower(opt$species) %in% c("human", "mouse") & is.null(opt$customref) ){
        print("特殊物种，未提供细胞类型鉴定数据集，不进行celltyping!")
    } else {
        if ( !is.null(opt$customref)){
            species = opt$species
            customref = opt$customref
        } else if ( tolower(opt$species) %in% c("human", "mouse") & is.null(opt$customref) ){
            species = opt$species
            print("人或小鼠未指定数据集，默认使用hpca或immgen数据集!")
            customref = ifelse(tolower(opt$species) == "human","/data/database/celltype_refdata/logNorm_rds/hpca.rds","/data/database/celltype_refdata/logNorm_rds/immgen.rds")
        }
        param = paste0("-i ",sub_h5seurat," -f h5seurat",
                       " -o ",root_dir,"/Reference_celltype  -d h5seurat  --assay ",opt$assay, 
                       " --dataslot counts --update T ",predicate," celltyping -r ",customref," --annolevel ", opt$annolevel,
                       " --usecluster F --demethod classic --pointsize ",opt$pointsize," -n 25 --reduct ",opt$reduct2," --species ",species)
        cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", param)
        cmd = c(cmd, "# 细胞类型鉴定", cmd1)
    }

}
## 基因可视化的命令组合  
if ( "marker_vis" %in% method ){
    target = basename(list.dirs("/data/database/sc_subtype_refmarker",recursive = FALSE, full.names = TRUE) )
    all_genelist = NULL
    if( is.null(opt$celltype) & is.null(opt$predicate) & is.null(opt$subset) & is.null(opt$extraGene)){
         print("既未指定细胞类型，又无法获取有效的细胞类型信息，不执行基因可视化。")
         system(paste0("touch ",root_dir,"/既未指定细胞类型，又无法获取有效的细胞类型信息，不执行基因可视化.txt"))
    } else {
        all_celltype = NULL
        if ( !is.null(opt$celltype)){
            all_celltype = unlist(strsplit(opt$celltype, ","))
        # 匹配筛选的细胞类型
        } else if ( !is.null(opt$predicate) | !is.null(opt$subset) ){
            all_celltype = all_types[!grepl("^\\d+$", all_types)] 
        }
        # 非常规物种不适用默认基因列表
        if (as.logical(opt$default_markers) & tolower(opt$species) %in% c("human", "mouse","rat")){
            wl = read.delim("/data/database/sc_subtype_refmarker/whitelist.xls",sep = "\t")
            whitelist = as.list(wl)
            for (i in all_celltype){
                # ## 模糊匹配
                # match_type = target[ agrep(i, target, max.distance = 0.2, ignore.case = TRUE)]
                # 白名单匹配
                matches = lapply(whitelist, function(x) i %in% x)
                match_type = names(which(unlist(matches)))
                # # 匹配到多个细胞类型，可能是T细胞和B细胞
                # if ( length(match_type) > 1 & all(c("T_cell", "B_cell") %in% match_type) ) {
                    # print(paste0("匹配到多个细胞类型：",paste(match_type,collapse = "、"),",进行二次匹配！"))
                    # if( grepl("t", tolower(i)) ){ match_type = "T_cell" } else if ( grepl("b", tolower(i)) ){ match_type = "B_cell" }
                    # all_genelist = c(all_genelist, list.files(paste0("/data/database/sc_subtype_refmarker/",match_type),recursive = FALSE, full.names = TRUE))
                # # 未匹配到细胞类型
                # } else 
                if ( length(match_type) == 0 ){ 
                    print(paste0("细胞类型：",match_type,"未匹配到对应文件！"))
                # 匹配到其他细胞类型
                } else {
                    print(paste0("匹配到细胞类型：",match_type))
                    genelists = list.files(paste0("/data/database/sc_subtype_refmarker/",match_type),recursive = FALSE, full.names = TRUE)
                    #选择组织器官
                    #需要组织器官的细胞类型 c("Fibroblast","Epithelial","Microglia","Neuron","Myeloid")
                    if( !is.null(opt$tissue) ){
                        if (any(grepl(opt$tissue, genelists))) { 
                            genelists = genelists[grep(opt$tissue, genelists)]
                        } else { 
                            print(paste0("未匹配到对应组织类型：",opt$tissue,"，未出基因可视化结果，请核查！"))
                            system(paste0("touch ",root_dir,"/未匹配到对应组织类型：",opt$tissue,"，未出基因可视化结果，请核查.txt"))
                            genelists = NULL 
                        }
                    } else {
                        print(paste0("未提供组织类型信息，提供非组织类型marker基因可视化，若无请核查。"))
                        system(paste0("touch ",root_dir,"/未提供组织类型信息，提供非组织类型marker基因可视化，若无请核查.txt"))
                        # 匹配非 reference_marker_数字_字符串.xls 格式的文件
                        tissue_pattern = "reference_marker_\\d+_[[:alpha:]]+\\.xls$"
                        genelists = genelists[!grepl(tissue_pattern, genelists)]
                    }
                    #选择物种
                    if( !is.null(opt$species) & tolower(opt$species) %in% c("human", "mouse") ){
                        if (any(grepl(opt$species, genelists))) { genelists = genelists[grep(opt$species, genelists)]} else { genelists = genelists }
                    }
                    all_genelist = c(all_genelist, genelists)
                }
            }
        } else { print("非常规物种，或者已指定不使用默认基因列表，可以自行指定基因列表进行可视化。") }
        ## 外部基因列表
        if ( !is.null(opt$extraGene)){
            extraGenes = unlist(strsplit(opt$extraGene, ","))
            all_genelist = c(all_genelist, normalizePath(extraGenes))
        }
        for (genes in all_genelist){
            basedir_name = unlist(strsplit(basename(dirname(genes)), "/"))[1]
            basefile_name = gsub("\\.[[:alnum:]]+$", "", basename(genes))
            genelist_name <- paste(basedir_name, basefile_name, sep="_")
            param = paste0("-i ",sub_h5seurat," -f h5seurat",
                           " -o ",root_dir,"/featureplot_vlnplot/",genelist_name,"  -d h5seurat  --assay ",opt$assay, 
                           " --dataslot data ",predicate," visualize -x ",genes," -g clusters --pointsize ",opt$pointsize,
                           " -m vlnplot,featureplot,dotplot --vcolors customecol2 --ccolors spectral --reduct ",opt$reduct2,
                           " --dotsplit T --dodge F ")
            cmd1 = paste0("Rscript /public/scRNA_works/pipeline/oesinglecell3/exec/sctool ", param)
            cmd = c(cmd, "# 基因可视化", cmd1)
        }
    
    }
}
## 运行命令
print("即将开始运行亚型分析，运行代码如下：")
cat(cmd, sep = "\n" )
writeLines(cmd,  file.path(root_dir,paste0(filename,"_subclusters.sh")))
system(paste0("bash ", file.path(root_dir,paste0(filename,"_subclusters.sh"))))

print("亚型分析完成！")
## 结果整理 ##
print("开始整理文件……")

# 整理clustering
if ( file.exists(paste0(root_dir,"/Clustering")) ){
    system(paste0("mkdir -p ",root_dir,"/",filename,"/Clustering"))
    system(paste0("cp -r ", root_dir,"/Clustering/",opt$reduct2,"_Dimension_Reduction/*  ", root_dir,"/",filename, "/Clustering"))
    system(paste0("cp -r ", root_dir,"/Clustering/*cluster*  ", root_dir,"/",filename, "/Clustering"))

    if ( opt$reduct1 != "pca"){
        system(paste0("cp -r ", root_dir,"/Clustering/sampleid-batchid.xls  ", root_dir,"/",filename, "/Clustering/"))
    }
}
# 整理marker
if ( file.exists(paste0(root_dir,"/Marker")) ){
    system(paste0("mkdir -p ",root_dir,"/",filename,"/Marker/top10_markers_visualize_for_each_cluster"))
    system(paste0("mv ",root_dir,"/Marker/markers_vis*  ",root_dir,"/",filename,"/Marker/top10_markers_visualize_for_each_cluster/"))
    # 获取以"anno.xls"结尾的文件列表
    anno_files <- list.files(path = paste0(root_dir,"/Marker/"), pattern = "anno.xls$", full.names = TRUE)
    if ( length(anno_files) > 0 ){
        system(paste0("cp ", root_dir,"/Marker/*_anno.xls  ", root_dir,"/",filename, "/Marker/"))
    } else {
        system(paste0("cp ", root_dir,"/Marker/*.xls  ", root_dir,"/",filename, "/Marker/"))
    }
    system(paste0("cp ", root_dir,"/Marker/topmarker_gene_heatmap*  ", root_dir,"/",filename, "/Marker/"))
}
# 整理细胞类型鉴定
if(file.exists(paste0(root_dir,"/Reference_celltype"))){
    system(paste0("mv ",root_dir,"/Reference_celltype  ",root_dir,"/",filename,"/"))
}
# 整理基因可视化
if(file.exists(paste0(root_dir,"/featureplot_vlnplot"))){
    system(paste0("mv ",root_dir,"/featureplot_vlnplot  ",root_dir,"/",filename,"/"))
}

print("搞定！")
