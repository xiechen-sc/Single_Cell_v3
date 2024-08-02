#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
# import
import os, re, sys
from single_cell_auto import *
# 获取命令行参数
args = sys.argv[1:]
config = args[0]
config_path = os.path.dirname(os.path.abspath(config))  # 生成文件与 config 文件同目录
# 参数
pattern = r"D.*OE.*\d+|HT\d+.*|ZO.*"
project_id = re.search(pattern,config_path)
if not project_id:
    project_id = "No Match"
else:
    project_id = project_id.group()
    project_id = re.sub('/.*','',project_id)

### yaml 转化为类属性
yaml_data = read_yaml_file(config)
# 获取具体分析的内容
module_analysis = yaml_data['run']
# 实例化并获取脚本
get_script_fun(module_analysis=module_analysis,config_path=config_path,yaml_data=yaml_data,project_id=project_id)





