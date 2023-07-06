from django.contrib import admin
from . import models
# Register your models here.
admin.site.register(models.unique_id)
admin.site.register(models.cookie_sequence)
admin.site.register(models.cookie_position)
admin.site.register(models.cookie_opedvar)
admin.site.register(models.cookie_sequence_index)
admin.site.register(models.cookie_position_index)
admin.site.register(models.cookie_opedvar_index)
