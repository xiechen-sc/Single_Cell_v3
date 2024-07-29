#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
# import
import yaml
import os, re, sys
import single_cell_auto
# 获取命令行参数
args = sys.argv[1:]
config = args[0]
config_path = os.path.dirname(os.path.abspath(config))  # 生成文件与 config 文件同目录
# 参数
pattern = r"D.*OE.*\d+|HT\d+.*"
project_id = re.search(pattern,config_path)
if not project_id:
    project_id = "No Match"
else:
    project_id = project_id.group()
    project_id = re.sub('/.*','',project_id)

# 定义函数 
# 读取 YAML 文件
def read_yaml_file(file_path):
    with open(file_path, 'r',encoding='utf-8') as file:
        yaml_data = yaml.safe_load(file)
    return yaml_data    
# 定义一个父类 保存公共数据
class ProjectInfo():

    outdir = config_path
    project_id = project_id
    ###  定义类 将 yaml 递归创建至 类中
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            if isinstance(value, dict):
                setattr(self, key, ProjectInfo(**value))
            else:
                setattr(self, key, value)
    
    def insert_sql(self):  # 增
        pass
    
    def delete_sql(self):  # 删
        single_cell_auto.delete(self.project_id)

    def updata_sql(self):  # 改
        pass

    def select_sql(self):  # 查
        return single_cell_auto.sql_scelect(self.project_id)

### yaml 转化为类属性
yaml_data = read_yaml_file(config)
# 获取具体分析的内容
module_analysis = yaml_data['run']

if module_analysis == 'featureplot':
    from scClass.featureplot import mkobj
    sc_obj = mkobj(father=ProjectInfo,yaml_data=yaml_data)
    print(sc_obj/)
    sc_obj.get_script()




