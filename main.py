#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
import os

# 参数
path = os.getcwd()
project_id = 'proj1234'
# 定义一个父类 保存公共数据
class ProjectInfo():

    outdir = path
    project_id = project_id

    def __init__(self):
        pass
    
    def insert_sql(self):  # 增
        pass
    
    def delete_sql(self):  # 删
        delete(self.project_id)

    def updata_sql(self):  # 改
        pass

    def select_sql(self):  # 查
        return sql_scelect(self.project_id)


temp = ProjectInfo()

print(temp.outdir)
ret = temp.select_sql()
print(ret)
