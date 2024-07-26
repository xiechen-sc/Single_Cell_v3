#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
# 用于获取远程数据库的数据包
import pymysql
def sql_scelect(project_id):
    # 数据库连接配置
    config = {
        'user': 'reader',        # MySQL 用户名
        'password': 'oe1311',    # MySQL 密码
        'host': '10.101.48.44', # 远程 MySQL 服务器的 IP 地址或主机名
        'database': 'test_sql',# 连接的数据库名
        'cursorclass': pymysql.cursors.DictCursor
    }

    # project_id = 'proj1234'  # 要查询的项目号

    try:
        # 建立数据库连接
        connection = pymysql.connect(**config)
        
        # 手动创建游标
        cursor = connection.cursor()
        try:
            sql = "SELECT * FROM cfg_info WHERE project_id = %s;"
            cursor.execute(sql, (project_id,))
            result = cursor.fetchone()

        finally:
            # 确保游标被关闭
            cursor.close()

    finally:
        if connection:
            connection.close()
    return result
if __name__ == '__main__':
    ret = sql_scelect('proj11234')
    print(ret)
