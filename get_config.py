#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
from single_cell_auto import *
from config_def import *
import sys
import os
def main():
    config_out = os.getcwd()
    if len(sys.argv) == 1: # 没有参数传入
        analysis_modules = show_guide()  # 打印交互信息
        analysis_modules = analysis_modules.split(',')
    else: # 有参数传入
        analysis_modules = sys.argv[1:]
        analysis_modules = analysis_modules[0].split(',')
    
    
    for analysis_module in analysis_modules:
        # 检查内部元素是否只包含 数字和 逗号
        if not analysis_module.isdigit():
            show_help(analysis_module=analysis_module)
        if analysis_module == '4':
            get_featureplot(config_out=config_out)
        elif analysis_module == '1':
            get_modified_cell_type(config_out=config_out)
        elif analysis_module == '3':
            get_diff_enrich(config_out=config_out)
        elif analysis_module == '2':
            get_sub_clusters(config_out=config_out)
        elif analysis_module == '5':
            get_singleR(config_out=config_out)
        elif analysis_module == '6':
            get_enrichment(config_out=config_out)
        elif analysis_module == '7':
            get_scenic(config_out=config_out)
        elif analysis_module == '0':
            get_decontX(config_out=config_out)
        elif analysis_module == '8':
            get_monocle2(config_out=config_out)
        elif analysis_module == '9':
            get_addmodulescore(config_out=config_out)
        else:
            print(f'{analysis_module} :error! No or to be developed')


if __name__ == '__main__':
    main()
