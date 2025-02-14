from .base_class import BaseClass
from single_cell_auto.util import cell_name_normalization,get_species_info,jinggao,database_add
from single_cell_auto.cmd_module import cmd_addmodulescore
# diff
class Addmodulescore(BaseClass):

    def get_script(self):
            analysis_module = self.run
            outdir = self.outdir
            # input = self.input # 输入 seurat 对象文件
            # species = self.species  # 物种
            # genelist = self.genelist  # 输入genelist 列表
            # groupby = self.groupby  # 分组方式
            # splitby = self.splitby  # 拆分展示小提琴图
            # fsplitby = self.fsplitby  # 拆分展示umap图
            # pvalue = self.pvalue
            # # 下方内容选择性填写 建议默认 
            # reduct = self.reduct # 降维方式 
            # output = self.output  # 输出目录
            # assay = self.reduct  # RNA  SCT
            # dataslot = self.dataslot  # counts,data,scale.data
            # strict = self.strict  # 是否使用严格模式筛选gene，默认FALSE
            # assay = self.assay
            
            cmd = cmd_addmodulescore(self)
            out_script = f'{outdir}/cmd_addmodulescore.sh'
            with open(out_script,"w") as f:
                f.write(cmd)
            print(f"脚本 {out_script} 已生成")
            #### 物种信息保存至数据库
            db_update_bg = self.pjif  
            db_update_bg['species'] = self.species
            self.update_info_bag = db_update_bg
                    

