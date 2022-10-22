from django.contrib import admin

# Register your models here.
from pegsub.models import *
from pegsub.models import *

# 注册站点模型
# admin.site.register(preInfo)
# admin.site.register(preFileInfo)
# admin.site.register(pegRNAInfo)
# #correct
# admin.site.register(Final_install_info)
# admin.site.register(Final_correct_info)


#模型类注册站点
admin.site.register(cookie_sequence)
admin.site.register(cookie_sequence_index)
admin.site.register(cookie_position)
admin.site.register(cookie_position_index)
admin.site.register(cookie_databse)
admin.site.register(cookie_database_index)
# 四月份新注册站点模型
admin.site.register(cookie_sequence_April)
admin.site.register(cookie_sequence_April_index)
#四月份新注册站点模型
admin.site.register(cookie_position_April)
admin.site.register(cookie_position_April_index)
#四月份新注册站点模型
admin.site.register(cookie_database_April)
admin.site.register(cookie_database_April_index)