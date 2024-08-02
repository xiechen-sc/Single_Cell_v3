#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
from single_cell_auto import *
from config_def import *
def main():
    import os
    config_out = os.getcwd()
    show_guide()  # 打印交互信息

    analysis_module = input("请选择分析内容!\n")
    # 判断输入是否正确
    # analysis_module = '4'


    if analysis_module == '4':
        get_featureplot(config_out)
    else:
        print('error! No or to be developed')
        exit(1)

if __name__ == '__main__':
    main()
