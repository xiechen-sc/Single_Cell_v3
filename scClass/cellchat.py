from .base_class import BaseClass
from single_cell_auto.util import cell_name_normalization,get_species_info,jinggao,database_add
# Cellchat
class Cellchat(BaseClass):

    def get_script(self):
            analysis_module = self.run
            version = self.version
            output = self.output
            input = self.input # 输入 seurat 对象文件
            species = self.species  # 物种
            column4cell = self.column4cell  # 用于进行通讯分析的metadata列名
            groupby = self.groupby  # 分组列名
            contrast = self.contrast  # 分组信息
            subsetby = self.subsetby  # 需要截取的列名
            which_cells = self.which_cells  # 需要截取的具体变量
            rds = self.rds  # 用于进行通讯分析的rds文件
            topn = self.topn  # 绘制气泡图筛选的top受配体对数
            strict = self.strict  # 是否使用严格模式筛选gene
            palette = self.palette  # 配色
            cellchatdb = self.cellchatdb  # cellchatdb
            # input 处理
            add_cmd,input = self.judgment_seurat_rds(input) # 这里返回 转化命令与新的 rds 路径
            # 简单处理物种
            if 'human' in species:
                species = 'human'
            elif 'mouse' in species:
                species = 'mouse'
            else:
            # 一般只能做人或者小鼠
            # if species not in ['mouse','human']:
                jinggao('特殊物种请注意！！！！！！！！！')
                jinggao('请 blast 或 手动制作数据库，并修改config 的物种参数！！！')
            if version == 'V2':
                # 处理 cmd 命令
                cmd = '### cellchat V2 ###\n'
                cmd += 'set -e\n'
                cmd += f'module purge && module load OESingleCell/v_3.0.0_visium_produce\n'
                cmd += add_cmd
                cmd += self.add_cmd_row(f'sctool cellchat2')
                cmd += self.add_cmd_row(f'--input {input}')
                cmd += self.add_cmd_row(f'--species {species}')
                cmd += self.add_cmd_row(f'--column4cell {column4cell}')
                cmd += self.add_cmd_row(f'--assay RNA')
                cmd += self.add_cmd_row(f'--palette {palette}')
                cmd += self.add_cmd_row(f'--topn {topn}')

                #  是否取子集，加上subsetby 和 which_cells
                if subsetby != "None":
                    cmd += self.add_cmd_row(f'--subsetby {subsetby}')
                    cmd += self.add_cmd_row(f'--which_cells {which_cells}')

                # 分组比较时，加上 groupby 和 contrast
                if groupby != "None":
                    cmd += self.add_cmd_row(f'--groupby {groupby}')
                    cmd += self.add_cmd_row(f'--contrast {contrast}')
                if strict != "F":
                    cmd += self.add_cmd_row(f'--strict {strict}')
                if rds != 'None':
                    cmd += self.add_cmd_row(f'--rds {rds}')
                cmd += self.add_cmd_row(f'--output {output}', end=True)
                
            elif version == 'V1':
                # 处理 cmd 命令
                cmd = '### cellchat V1 ###\n'
                cmd += 'set -e\n'
                cmd += f'module purge && module load OESingleCell/3.0.d\n'
                cmd += add_cmd
                cmd += self.add_cmd_row(f'Rscript /public/scRNA_works/pipeline/scRNA-seq_further_analysis/CellChat_v1.6.1.R')
                cmd += self.add_cmd_row(f'-i {input}')
                cmd += self.add_cmd_row(f'-s {species}')
                cmd += self.add_cmd_row(f'-c {column4cell}')
                cmd += self.add_cmd_row(f'--assay RNA')
                cmd += self.add_cmd_row(f'--palette {palette}')
                cmd += self.add_cmd_row(f'--topn {topn}')

                #  是否取子集，加上subsetby 和 which_cells
                if subsetby != "None":
                    cmd += self.add_cmd_row(f'-q {subsetby}')
                    cmd += self.add_cmd_row(f'-u {which_cells}')

                # 分组比较时，加上 groupby 和 contrast
                if groupby != "None":
                    cmd += self.add_cmd_row(f'-g {groupby}')
                    cmd += self.add_cmd_row(f'-d {contrast}')
                if strict != "F":
                    cmd += self.add_cmd_row(f'--strict {strict}')
                if rds != 'None':
                    cmd += self.add_cmd_row(f'--rds {rds}')
                
                if cellchatdb != 'None':
                    cmd += self.add_cmd_row(f'--cellchatdb {cellchatdb}')
                cmd += self.add_cmd_row(f' -o {output}', end=True)
            else:
                print("版本错误，请检查！")

            # out_script="run_cellchat.sh"
            out_script = f'{self.outdir}/cmd_{self.run}.sh'
            with open(out_script, "w") as f:
                f.write(cmd)
            print(f"脚本{out_script}已生成")
