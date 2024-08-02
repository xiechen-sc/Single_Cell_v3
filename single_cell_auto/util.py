puple = '\033[35m'
reset = '\033[0m'
cyan = '\033[36m'
yellow = '\033[33m'
red = '\033[31m'
# 交互信息
def show_guide():  
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
    ipt = input(f"{yellow}请选择分析内容!一次选多个请使用 ',' 分隔！{reset}\n : ")
    return ipt
    
# 读取 YAML 文件
def read_yaml_file(file_path):
    import yaml
    with open(file_path, 'r',encoding='utf-8') as file:
        yaml_data = yaml.safe_load(file)
    return yaml_data 

def show_help(analysis_module=None):
    print('###################################################################################################################################')
    print(f"输入 {analysis_module} 错误 : {cyan}输入的字符串，只能是 数字 和 逗号组成！数字表示分析类型，多个分析用逗号分隔！！！{reset}")  
    print('###################################################################################################################################')
    exit(1)
    
# 获取 物种 配置信息
def get_species_info(species=None):
    species_info_file = '/gpfs/oe-scrna/zhengfuxing/Single_Cell_v3/config_def/species_config.yaml'
    species_dict = read_yaml_file(species_info_file)
    try:
        species_info = species_dict[species]
    except KeyError:
        return 'None'
    else:
        return species_info
    
# 警告信息
def jinggao(str):
    info = f"{red}{str}{reset}"
    print(info)







if __name__ == '__main__':
    info = get_species_info('asd')
    print(info)