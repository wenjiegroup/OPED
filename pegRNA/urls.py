from pegRNA import views
from django.urls import path

urlpatterns = [
    path("", views.index, name="index"),

    path("Sequence/", views.page_sequence, name="sequence"),
    path("Position/", views.page_position, name="position"),
    path("OPEDVar/", views.page_opedvar, name="opedvar"),

    path("Sequence/history/", views.history_page, name="seq_his"),
    path("Position/history/", views.history_page, name="pos_his"),
    path("OPEDVar/history/", views.history_page, name="opedvar_his"),

    path("Sequence/download/", views.download_page, name="seq_dl"),
    path("Position/download/", views.download_page, name="pos_dl"),
    path("OPEDVar/download/", views.download_page, name="opedvar_dl"),

    path("Sequence/history/detail/", views.detail_page, name="detail_seq"),
    path("Position/history/detail/", views.detail_page, name="detail_pos"),
    path("OPEDVar/history/detail/", views.detail_page, name="detail_opedvar"),

    path("help/", views.page_help, name="help"),
    path("about/", views.page_about, name="about"),

]
