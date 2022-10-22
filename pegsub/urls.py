from django.conf.urls import url
from django.http import JsonResponse
from django.shortcuts import render

from pegsub import views,djangohtml
# from pegsub.models import pegRNAInfo

from django.urls import include, path

urlpatterns = [
    path("", views.index, name="index"),
    path("home/", views.to_home, name="home"),
    path("pegRNA_PredictingCodes/", views.add_predict_show, name="qsub_Data"),

    #sequence

    path("pegRNA_file/", views.add_file_predict_show, name="qsub_file"),

    path("downloadfile_update/", views.downloadfile_update,name="downloadfile_update"),
    path("downloadfile_update_history/", views.downloadfile_update_history, name="downloadfile_update_history"),

    path("pegRNA_file_update/",views.add_file_predict_show_update,name="qsub_file_update"),

    #pos
    path("pegRNA_file_pos/",views.add_file_predict_show_pos,name="qsub_file_pos"),

    path("downloadfile_position/",views.downloadfile_position,name="downloadfile_position"),
    path("downloadfile_position_history/", views.downloadfile_position_history, name="downloadfile_position_history"),



    # '''
    # 四月更新为OPEDVar
    # '''
    path("OPEDVar/",views.search_database,name="OPEDVar"),
    #path("OPEDVar/", views.OPEDVar, name="OPEDVar"),

    path("downloadfile_database/", views.downloadfile_database, name="downloadfile_database"),
    path("downloadfile_database_history/", views.downloadfile_database_history, name="downloadfile_database_history"),



    #shell from liufeng


    # path("Simple_Version/",views.Simple_Version,name="simple"),

##ajax_try

    # path("ajax_predict/",views.ajax_input,name="ajax_input"),
    path("to_predict/",views.to_predict,name="to_predict"),
    path("predict_ajax/",views.ajax_predict),

    # path("history/",views.request_update),

    path("history_update_index/",views.request_update,name="requsert_update_index"),
    path("detail_sequence/", views.detail_sequence, name="detail_sequence"),


    path("history_position_index/", views.request_position, name="requsert_position_index"),
    path("detail_position/", views.detail_position, name="detail_sequence"),
    #=Simple_version

    path("history_database_index/", views.request_database, name="requsert_position_index"),
    path("detail_database/", views.detail_database, name="detail_sequence"),
    #分页
    # path("peg_list/", views.pegRNA_list, name="peg_list"),
    # 查询分页路由
    # 主页面路由


]
