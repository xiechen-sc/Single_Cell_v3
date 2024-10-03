import os
import yaml
import re
puple = '\033[35m'
reset = '\033[0m'
cyan = '\033[36m'
yellow = '\033[33m'
red = '\033[31m'
# 处理特殊字符
def cell_name_normalization(cell_name):
    pattern = r'\s+|\/|\(|\)|\||\&|\^|\%|,'
    new_name = re.sub(pattern, '_', cell_name)
    return new_name
# 判断项目号
def get_project_id(config_path):
    pattern = r"D[^J].*OE.*\d+|HT\d+.*|ZO.*|DZQD.*"
    project_id = re.search(pattern,config_path)
    if not project_id:
        project_id = "No Match"
    else:
        project_id = project_id.group()
        project_id = re.sub('/.*','',project_id)
    return project_id
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
        5.singleR 注释细胞类型
        6.enrichment 富集分析
        {reset}
        """
            )
    ipt = input(f"{yellow}请选择分析内容!一次选多个请使用 ',' 分隔！{reset}\n : ")
    return ipt
    
# 读取 YAML 文件
def read_yaml_file(file_path):
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

# 分布式数据库
# 获取数据库创建地址
def get_database_path(config_path,project_id):
    if project_id == "No Match":
        jinggao("获取project_id 变量时失败 无法创建或检索数据库，请手动填写 config ！！！")
        database_path = 'None'
        return 
    else:
        database_path = config_path.split(project_id)[0] + project_id + '/.project_info'
        return database_path
# 写入数据库文件
def save_dict_to_yaml(data_base_file, project_info):
        try:
            with open(data_base_file, 'w') as file:
                yaml.dump(project_info, file, default_flow_style=False, allow_unicode=True)
        except PermissionError:
            print(f"该项目的根目录下没有权限写入，请联系所有人开放写入权限！否则影响自动化工作！")

# 查询 一般发生在 已有数据库的情况下 去快速获取数据
def database_retrieval(config_path):
    project_id = get_project_id(config_path)
    data_base_file = get_database_path(config_path=config_path,project_id=project_id)
    if os.path.exists(data_base_file):
        project_info = read_yaml_file(data_base_file)
        return project_info
    else:
        return {} # 保证查询返回结果一定是字典

# 创建 workflow 之后执行 选择需要上传数据库的部分执行
def database_add(config_path,config_info):  # 传入的 config_info 是一个字典  由模块开发者认为有必要添加进数据库的内容 会往其中添加 
    project_id = get_project_id(config_path)
    data_base_file = get_database_path(config_path=config_path,project_id=project_id)
    if os.path.exists(data_base_file):  # 数据库存在
        project_info = read_yaml_file(data_base_file)
        project_info = project_info | config_info
        
    else: # 不存在数据库 
        project_info = {}
        project_info['project_id'] = project_id
        project_info = project_info | config_info

    save_dict_to_yaml(data_base_file=data_base_file, project_info=project_info)
    
def v4tov3():
    pass # 待良钰开发

def Frequent_species():  # 记录常见的物种
    fs = ['human','human_2020','human_2024','mouse','mouse_2020','mouse_2024']
    return(fs)




if __name__ == '__main__':
    info = get_species_info('asd')
    print(info)