#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
def main():
    import os
    #  color
    puple = '\033[35m'
    reset = '\033[0m'
    cyan = '\033[36m'
    yellow = '\033[33m'
    # 得到当前路径
    config_out = os.getcwd()

    print(f"""{cyan}
    ################################################################################################################################
    #########################                    欢迎使用单细胞自动化售后工具                      #################################
    ################################################################################################################################
    {reset}
    {yellow}
    1.修改细胞类型
    2.亚群分析
    3.差异分析
    4.marker基因绘图



    {reset}
    """
        )

    analysis_module = input("请选择分析内容!\n")
    # 判断输入是否正确
    # analysis_module = '4'


    if analysis_module == '4':
        from config_def import marke_featrueplot
        marke_featrueplot.get_featureplot(config_out)
    else:
        print('error! No or to be developed')
        exit(1)

if __name__ == '__main__':
    main()
