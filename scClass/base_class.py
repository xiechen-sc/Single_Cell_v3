from single_cell_auto import *
# 定义一个父类 保存公共数据
class BaseClass:

    ###  定义类 将 yaml 递归创建至 类中
    def __init__(self,config_path=None,project_id=None,**kwargs):
        self.outdir = config_path
        self.project_info = self.select_sql()
        self.project_id = project_id
        for key, value in kwargs.items():
            if isinstance(value, dict):
                setattr(self, key, BaseClass(**value))
            else:
                setattr(self, key, value)

    
    def insert_sql(self):  # 增
        pass
    
    def delete_sql(self):  # 删
        delete(self.project_id)

    def updata_sql(self):  # 改
        pass

    def select_sql(self):  # 查
        project_info = database_retrieval(self.outdir)
        return project_info
    
    def add_cmd_row(self,cmd,end=False):
        if end:
            normalize_cmd = cmd + '\n\n'
        else:
            normalize_cmd = cmd + ' \\\n'
        return normalize_cmd
