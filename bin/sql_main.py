#!/gpfs/oe-scrna/zhengfuxing/conda/scAuto/bin/python
from Select_sql import sql_scelect

ret = sql_scelect('proj123')
if not ret:
    print('var is not ')
else:
    print(ret)

