from single_cell_auto import *
import re,ast
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
        self.get_project_info()
        self.update_info_bag = {}

    # 程序结束更新数据库
    def __del__(self):
        self.project_info_update()

    
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
    
    def seurat2rds(self,seurat,outdir):
        outdir2 = outdir + '/rds'
        outdir3 = outdir2  # 初始化 outdir3
        counts = 0
        while os.path.exists(outdir3):
            outdir3 = outdir2 + f'_v{counts}'
            counts += 1
        outdir2 = outdir3
        out_rds = outdir2 + '/data_ob_v3.rds'
        cmd_seurat2rds = f"\n########## seurat2rds\n/gpfs/oe-scrna/guopengyu/script/rds.sh -i {seurat} -o {outdir2}\n########## seurat2rds finish!\n\n"
        
        return cmd_seurat2rds,out_rds
    
    def judgment_seurat_rds(self,input_file):
        outdir = self.outdir
        lst = input_file.split('.')
        lst_len = len(lst)
        suffix = lst[lst_len-1]
        if suffix == 'h5seurat':
            cmd_seurat2rds,out_rds = self.seurat2rds(seurat=input_file,outdir=outdir)
            return cmd_seurat2rds,out_rds
        elif suffix == 'rds':
            return '',input_file
        else:
            jinggao('输入文件后缀错误，只能是 rds 或 h5seurat!')
            exit()
    
    # 数据库信息添加至对象的属性
    def get_project_info(self):
        config_path = self.outdir
        self.pjif = database_retrieval(config_path)

    # 更新数据库信息
    def project_info_update(self):  # 传入的 update_info_bag 是一个字典  由模块开发者认为有必要添加进数据库的内容 会往其中添加 
        self.pjif = self.pjif | self.update_info_bag
        config_path = self.outdir
        project_id = get_project_id(config_path)
        data_base_file = get_database_path(config_path=config_path,project_id=project_id)
        # print(self.pjif)
        save_dict_to_yaml(data_base_file=data_base_file, project_info=self.pjif)

    # 为输入的字符串列表添加引号！ 并将字符串转化为列表后返回
    def quotation_mark(self,raw_str):
        mid_str = raw_str.replace('\"','').replace("\'",'')
        new_str = re.sub(r'([^,\[\]]+)',r'"\1"',mid_str)
        lst = ast.literal_eval(new_str)
        lst2 = self.get_select(lst)
        return lst2


    def get_select(self,ipt):
        if isinstance(ipt, list):
            lst2 = []
            for i in ipt:
                ret = self.get_select(i)
                lst2.append(ret)
            return lst2
        else:
            ret = re.sub(r'\s+$|^\s+','',ipt) 
            return(ret)
        
