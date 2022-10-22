#!/usr/bin/python3
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, FileResponse, StreamingHttpResponse, JsonResponse, request

from pegRNA_PredictingCodes.predict_efficiency_of_ClinVar import read_data_of_Single
from pegsub import models, hjs_change
from pegRNA_PredictingCodes import predict_efficiency_of_ClinVar, train_model, evaluate_model
import torch
import pandas as pd
import json
from . import views
from django.contrib import messages
import datetime
from decimal import Decimal
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

# 分页用到的import
# from pegsub.models import pegRNAInfo
from django.core.paginator import Paginator, PageNotAnInteger, InvalidPage, EmptyPage


# 导入预测的py文件

# Create your views here.

# 定向网页
def index(request):
    return render(request, "index.html")


# 定向简单版本
def Simple_Version(request):
    return render(request, "Simple_version.html")


# deredict index
def to_home(request):
    return render(request, "index.html")


# 提交预测字段
def add_predict_show(request):
    if request.method == "POST":
        Target = request.POST.get("Target", None)
        PBS = request.POST.get("PBS", None)
        RT = request.POST.get("RT", None)
        result = hjs_change.hjs_test(Target, PBS, RT)
        # models.preInfo.objects.create(
        #     Target=Target,
        #     PBS=PBS,
        #     RT=RT,
        #     result=result,
        # )
        result = result.replace('[', '').replace(']', '')
        print(result)
        return render(request, "predict.html", {"result": result})
    return render(request, "predict.html")


# direct——predict
def to_predict(request):
    return render(request, "predict.html")


# AJAX
@csrf_exempt
def ajax_input(request):
    if request.method == "POST":
        response = HttpResponse()
        Target = request.POST.get("Target", None)
        PBS = request.POST.get("PBS", None)
        RT = request.POST.get("RT", None)
        print(Target, PBS, RT)
        result = hjs_change.hjs_test(Target, PBS, RT)
        return response


@csrf_exempt
def ajax_predict(request):
    #     if request.method == "POST":
    #         response=HttpResponse()
    #         Target=
    pass


# 提交文件进行预测
def add_file_predict_show(request):
    if request.method == "POST":
        seq = request.POST.get("seq_input", None)
        # print(type(seq))
        info_list = hjs_change.hjs_file(seq)
        # info_list=seq
        print(info_list)
        return render(request, "predict_file.html", {"info_list": info_list})

    return render(request, "predict_file.html")


# 让info列返回
# score小数点后两位
#（1）第一种
def decimal2Effi(x):
    return Decimal(x['Efficiency']).quantize(Decimal("0.00"))
#（2）第二种
def decimal2(x):
    return Decimal(x['score']).quantize(Decimal("0.00"))



# 颜色不一样的返回页面
def getTarget1(x):
    return x['Target'][:4]


def getSpacer(x):
    return x['Target'][4:24]


def getNGG(x):
    return x['Target'][24:27]


def getlast_NGG(x):
    return x['Target'][27:]


def getNG(x):
    return x['Target'][24:26]


def getlast_NG(x):
    return x['Target'][26:]


# qsub file to predict update
def add_file_predict_show_update(request):
    if request.method == "POST":
        user_seq = request.POST.get("seq_input", None)
        pam_type = request.POST.get("pam_type", None)
        dis_nickase = request.POST.get("dis_nickase", None)
        max_len_PBS = request.POST.get("max_len_PBS", None)
        max_len_RT = request.POST.get("max_len_RT", None)
        min_len_PBS = request.POST.get("min_len_PBS", None)
        min_len_RT = request.POST.get("min_len_RT", None)
        top_n = request.POST.get("top_n", None)
        #April 更新
        min_dis_sgRNA_nickase=request.POST.get("min_dis_sgRNA_nickase",None)
        max_dis_sgRNA_nickase=request.POST.get("max_dis_sgRNA_nickase",None)
        homology=request.POST.get("homology",None)
        # 获取提交计算的当前时间
        curr_time = datetime.datetime.now()
        time_str = curr_time.strftime('%Y%m%d%H%M%S')

        input_dict = {
            'user_seq': user_seq,
            'pam_type': pam_type,
            'dis_nickase': dis_nickase,
            'max_len_PBS': str(max_len_PBS),
            'max_len_RT': str(max_len_RT),
            'min_len_PBS': str(min_len_PBS),
            'min_len_RT': str(min_len_RT),
            'top_n': top_n,
            'min_dis_sgRNA_nickase':str(min_dis_sgRNA_nickase),
            'max_dis_sgRNA_nickase':str(max_dis_sgRNA_nickase),
            'homology':str(homology)

        }
        print(input_dict)

        info_list, ERROR = hjs_change.hjs_file_update(input_dict, time_str)

        info_list['score'] = info_list.apply(decimal2, axis=1)

        print("houtaiyunxing_update")
        print(info_list)
        """ 此注释部分为更新前update的Target 分颜色部分
        if ERROR:
            pass
        else:
            # 区分颜色
            info_list.loc[:, 'Target1'] = info_list.apply(getTarget1, axis=1)
            info_list.loc[:, 'Target2'] = info_list.apply(getSpacer, axis=1)
            print(pam_type)
            if pam_type == "NGG":
                info_list.loc[:, 'Target3'] = info_list.apply(getNGG, axis=1)
                info_list.loc[:, 'Target4'] = info_list.apply(getlast_NGG, axis=1)
            else:
                info_list.loc[:, 'Target3'] = info_list.apply(getNG, axis=1)
                info_list.loc[:, 'Target4'] = info_list.apply(getlast_NG, axis=1)

        print(info_list)
        """

        print("ER",ERROR)
        start = 1
        return render(request, "predict_file_update.html", {"info_list": info_list, "ERROR": ERROR, "start": start})
    return render(request, "predict_file_update.html")


###qsub pos_file
def add_file_predict_show_pos(request):
    if request.method == "POST":
        chr = request.POST.get("chr", None)
        start_position = request.POST.get("start_position", None)
        edit_pattern = request.POST.get("edit_pattern", None)
        pam_type = request.POST.get("pam_type", None)
        dis_nickase = request.POST.get("dis_nickase", None)
        max_len_PBS = request.POST.get("max_len_PBS", None)
        max_len_RT = request.POST.get("max_len_RT", None)
        min_len_PBS = request.POST.get("min_len_PBS", None)
        min_len_RT = request.POST.get("min_len_RT", None)
        top_n = request.POST.get("top_n", None)

        #April 更新
        min_dis_sgRNA_nickase=request.POST.get("min_dis_sgRNA_nickase",None)
        max_dis_sgRNA_nickase=request.POST.get("max_dis_sgRNA_nickase",None)
        homology=request.POST.get("homology",None)

        ###
        # 获取系统时间即提交数据时间
        curr_time = datetime.datetime.now()
        time_str = curr_time.strftime('%Y%m%d%H%M%S')

        input_dict = {

            'chr': chr,
            'start_position': start_position,
            'edit_pattern': edit_pattern,
            'pam_type': pam_type,
            'dis_nickase': dis_nickase,
            'min_dis_sgRNA_nickase': str(min_dis_sgRNA_nickase),
            'max_dis_sgRNA_nickase': str(max_dis_sgRNA_nickase),
            'max_len_PBS': str(max_len_PBS),
            'max_len_RT': str(max_len_RT),
            'min_len_PBS': str(min_len_PBS),
            'min_len_RT': str(min_len_RT),
            'homology':str(homology),
            'top_n': top_n,


        }
        print(input_dict)
        # 有可能出现表没有结果的情况
        # Whether or not an error occurs is encapsulated in the function of hjs_file_pos
        info_list, ERROR = hjs_change.hjs_file_pos(input_dict, time_str)
        if ERROR=="None ERROR":

            info_list['score'] = info_list.apply(decimal2, axis=1)

        else:
            pass

        '''不显示targrt 此注释部分update之前Target部分
        if ERROR:
            pass
        else:
            # 区分颜色
            info_list.loc[:, 'Target1'] = info_list.apply(getTarget1, axis=1)
            info_list.loc[:, 'Target2'] = info_list.apply(getSpacer, axis=1)
            print(pam_type)
            if pam_type == "NGG":
                info_list.loc[:, 'Target3'] = info_list.apply(getNGG, axis=1)
                info_list.loc[:, 'Target4'] = info_list.apply(getlast_NGG, axis=1)
            else:
                info_list.loc[:, 'Target3'] = info_list.apply(getNG, axis=1)
                info_list.loc[:, 'Target4'] = info_list.apply(getlast_NG, axis=1)
        '''
        print("houtaiyunxing")


        start = 1
        return render(request, "predict_file_pos.html", {"info_list": info_list, 'ERROR': ERROR, 'start': start})
    return render(request, "predict_file_pos.html")


# search_Data
def search_database(request):
    if request.method == "POST":
        #四月更新字段
        queryType=request.POST.get("query_type",None)
        #四月更新字段end
        # alleleID = request.POST.get("query", None)
        queryItem=request.POST.get("queryItem",None)
        pam_type = request.POST.get("pam_type", None)
        direction = request.POST.get("direction", None)

        input_dict = {
            # 四月更新字段
            'queryType':queryType,
            # 四月更新字段end
            'queryItem': str(queryItem),
            'pam_type': pam_type,
            'direction': direction
        }

        curr_time = datetime.datetime.now()
        time_str = curr_time.strftime('%Y%m%d%H%M%S')

        print(input_dict)
        info_list, ERROR = hjs_change.hjs_file_database(input_dict, time_str)

        print(info_list)

        start = 1
        return render(request, "data_serach.html", {"info_list": info_list, "ERROR": ERROR, "start": start})
    return render(request, "data_serach.html")


# 文件下载###################################

def downloadfile_update(request):
    # file = open(f'Sequence/show_pegRNA_update.User.csv', 'rb')

    #四月份更新
    file=open(f'Sequence/show_pegRNA_update_April.csv','rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_sequence_result.csv"'
    return response


def downloadfile_update_history(request):
    file = open(f'Sequence/show_pegRNA_sequence_April_history_temp.User.csv','rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_sequence_result_history.csv"'
    return response


def downloadfile_position(request):
    # file = open(f'Sequence/show_pegRNA_position.User.csv', 'rb')
    file=open(f'Position/show_pegRNA_position_April.User.csv','rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_position_result.csv"'
    return response


def downloadfile_position_history(request):
    # file = open(f'Sequence/show_pegRNA_position.User.csv', 'rb')
    file=open(f'Position/show_pegRNA_position_April_history_temp.User.csv','rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_position_result.csv"'
    return response




def downloadfile_database(request):
    # file = open(f'Data/show_pegRNA_database.User.csv', 'rb')
    file = open(f'Data/show_pegRNA_database_April.User.csv', 'rb')
    response = FileResponse(file)
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_database_result.csv"'
    return response


def downloadfile_database_history(request):
    # file = open(f'Data/show_pegRNA_database_history_temp.User.csv')
    file = open(f'Data/show_pegRNA_database_history_April_temp.User.csv','rb')
    response = FileResponse(file)
    response['Content-Type'.encode('utf')] = 'application/octet-stream'
    response['Content-Disposition'] = 'attachment;filename="pegRNA_database_result_history.csv"'
    return response


# 历史记录###############################################

# # def history_update(request):
#     if request.method == "POST":
#         all = models.cookie_sequence.objects.all().order_by("-Time")
#         start_history = 1
#         return render(request, "predict_file_update.html", {"histroy_update_df": all, "start_history": start_history})
#     return render(request, "predict_file_update.html")


def request_update(request):
    # 实现搜索
    # key = request.GET.get('key')

    # all_users = models.cookie_sequence_index.objects.all().order_by("-time")
    #更新为四月份
    all_users=models.cookie_sequence_April_index.objects.all().order_by("-time")
    if all_users:
        paginator = Paginator(all_users, 10)
        page = request.GET.get('page')
        try:
            contacts = paginator.page(page)
            print(contacts)
        except PageNotAnInteger:
            contacts = paginator.page(1)
        except EmptyPage:
            contacts = paginator.page(paginator.num_pages)

        start_history = 1
        print(contacts)
        return render(request, 'predict_file_update.html', {'contacts': contacts, 'start_history': start_history})

    else:
        start_history = 1
        history_updateinfo = '暂无历史数据'
        return render(request, 'predict_file_update.html', {"history": history_updateinfo})


def detail_sequence(request):
    print(request.GET)
    time_id = request.GET.get("time_id", None)
    pam_type = request.GET.get("pam_type", None)
    '''此部分为更新前部分
    filter_model = models.cookie_sequence.objects.filter(time=time_id)
    filter_model.df = pd.DataFrame(list(models.cookie_sequence.objects.filter(time=time_id).values()))
    print(filter_model.df)
    filter_model.df.loc[:, 'Target1'] = filter_model.df.apply(getTarget1, axis=1)
    filter_model.df.loc[:, 'Target2'] = filter_model.df.apply(getSpacer, axis=1)
    filter_model.df['Efficiency']=filter_model.df.apply(decimal2Effi,axis=1)
    print(pam_type)
    if pam_type == "NGG":
        filter_model.df.loc[:, 'Target3'] = filter_model.df.apply(getNGG, axis=1)
        filter_model.df.loc[:, 'Target4'] = filter_model.df.apply(getlast_NGG, axis=1)
    else:
        filter_model.df.loc[:, 'Target3'] = filter_model.df.apply(getNG, axis=1)
        filter_model.df.loc[:, 'Target4'] = filter_model.df.apply(getlast_NG, axis=1)
    print(filter_model.df)
    '''

    '''此部分为四月更新后部分'''
    filter_model=models.cookie_sequence_April.objects.filter(time=time_id)
    filter_model.df = pd.DataFrame(list(models.cookie_sequence_April.objects.filter(time=time_id).values()))



    '''更新后部分end'''
    # filter_model.df.to_csv(f'Sequence/show_pegRNA_sequence_history_temp.User.csv')
    filter_model.df.to_csv(f'Sequence/show_pegRNA_sequence_April_history_temp.User.csv')
    start_history = 2
    # print(filter_model.df['Target3'])
    return render(request, 'predict_file_update.html', {'contacts_df': filter_model.df, 'start_history': start_history})
    return render(request, "predict_file_update.html")




def request_position(request):
    # 实现搜索
    # key = request.GET.get('key')

    all_users = models.cookie_position_April_index.objects.all().order_by("-time")

    if all_users:
        paginator = Paginator(all_users, 10)
        page = request.GET.get('page')
        try:
            contacts = paginator.page(page)
            print(contacts)
        except PageNotAnInteger:
            contacts = paginator.page(1)
        except EmptyPage:
            contacts = paginator.page(paginator.num_pages)

        start_history = 1
        print(contacts)
        return render(request, 'predict_file_pos.html', {'contacts': contacts, 'start_history': start_history})

    else:
        start_history = 1
        history_updateinfo = '暂无历史数据'
        return render(request, 'predict_file_pos.html', {"history": history_updateinfo})


def detail_position(request):
    print(request.GET)
    time_id = request.GET.get("time_id", None)
    pam_type = request.GET.get("pam_type", None)
    filter_model = models.cookie_position.objects.filter(time=time_id)
    filter_model.df=pd.DataFrame(list(models.cookie_position_April.objects.filter(time=time_id).values()))

    filter_model.df.to_csv(f'Position/show_pegRNA_position_April_history_temp.User.csv')
    '''此部分为四月更新前部分
    filter_model.df = pd.DataFrame(list(models.cookie_position.objects.filter(time=time_id).values()))
    print(filter_model.df)

    filter_model.df.loc[:, 'Target1'] = filter_model.df.apply(getTarget1, axis=1)
    filter_model.df.loc[:, 'Target2'] = filter_model.df.apply(getSpacer, axis=1)
    filter_model.df['Efficiency'] = filter_model.df.apply(decimal2Effi, axis=1)
    print(pam_type)
    if pam_type == "NGG":
        filter_model.df.loc[:, 'Target3'] = filter_model.df.apply(getNGG, axis=1)
        filter_model.df.loc[:, 'Target4'] = filter_model.df.apply(getlast_NGG, axis=1)
    else:
        filter_model.df.loc[:, 'Target3'] = filter_model.df.apply(getNG, axis=1)
        filter_model.df.loc[:, 'Target4'] = filter_model.df.apply(getlast_NG, axis=1)

    print(filter_model.df)
    filter_model.df.to_csv(f'Position/show_pegRNA_position_history_temp.User.csv')
    
    print(filter_model.df['Target3'])
    '''
    start_history = 2
    return render(request, 'predict_file_pos.html', {'contacts_df': filter_model.df, 'start_history': start_history})
    return render(request, "predict_file_pos.html")


def request_database(request):
    # 实现搜索
    # key = request.GET.get('key')

    all_users = models.cookie_database_April_index.objects.all().order_by("-time")

    if all_users:
        paginator = Paginator(all_users, 10)
        page = request.GET.get('page')
        try:
            contacts = paginator.page(page)
            print(contacts)
        except PageNotAnInteger:
            contacts = paginator.page(1)
        except EmptyPage:
            contacts = paginator.page(paginator.num_pages)

        start_history = 1
        print(contacts)
        return render(request, 'data_serach.html', {'contacts': contacts, 'start_history': start_history})

    else:
        start_history = 1
        history_updateinfo = '暂无历史数据'
        return render(request, 'data_serach.html', {"history": history_updateinfo})


def detail_database(request):
    print(request.GET)
    time_id = request.GET.get("time_id", None)

    # filter_model = models.cookie_databse.objects.filter(time=time_id)
    # filter_model.df = pd.DataFrame(list(models.cookie_databse.objects.filter(time=time_id).values()))
    '''四月更新'''
    filter_model = models.cookie_database_April.objects.filter(time=time_id)
    filter_model.df = pd.DataFrame(list(models.cookie_database_April.objects.filter(time=time_id).values()))

    # filter_model.df.to_csv(f'Sequence/show_pegRNA_database_history_temp.User.csv')
    filter_model.df.to_csv(f'Data/show_pegRNA_database_history_April_temp.User.csv')
    start_history = 2
    # print(filter_model.df['Target3'])
    return render(request, 'data_serach.html', {'contacts_df': filter_model.df, 'start_history': start_history})
    return render(request, "data_serach.html")


# def request_update_index(request):
#     #从URL中获取参数
#     page_num=request.Get.get("page")
#     print(page_num,type(page_num))
#     page_num=int(page_num)
#     #定义两个变量保存数据从哪取到哪
#     data_start = (page_num - 1) * 10
#     data_end = page_num * 10
#     #总数
#     all_users_count = models.cookie_sequence_index.objects.all().order_by("-time").count()
#     per_page = 10
#     #总共需要多少页码来显示
#     total_page,m=divmod(all_users_count,per_page)
#     # 页面上最多展示的页码
#     max_page = 11
#     half_max_page = max_page // 2
#     # 页面上展示的页码的开始页
#     page_start = page_num - half_max_page
#     # 页面上展示的页码的结束页
#     page_end = page_num + half_max_page
#     # 如果当前页减一半比 1 小
#     if page_start <= 1:
#         page_start = 1
#         page_end = max_page
#     # 如果当前页加一半比总页码还大
#     if page_end > total_page:
#         page_end = total_page
#         page_start = total_page - max_page + 1
#     if m:
#         total_page += 1
#     all_user=models.cookie_sequence_index.objects.all().order_by("-time")[data_start:data_end]
#
#
#
#
#     pass

# 提交字段进行预测
# def predict_ClinVar_single(gpu_id=1,type="Insertion",):
#     device = torch.device(f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu")
#     transformer = train_model.load_model(device, model_dir='Model_Trained',
#                                          model_name='pegRNA_Model_Merged_saved.order3_decoder_ori.pt')  # load model
#
#     o = evaluate_model.transformer_predictor_order3(transformer, data_X, 1024, device)[0]

def reverse_seq(seq):
    """get the reversed sequence of the input sequence
    Args:
        seq: the input sequence
    Returns:
        reversed sequence
    """

    return "".join(reversed(list(seq)))


def complement_seq(seq):
    """get the complementary sequence of the input sequence
    Args:
        seq: the input sequence
    Returns:
        complementary sequence
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([complement_map[s.upper()] for s in list(seq)])


# django数据库分页


###django数据库分页
# def info_list(request):
#     pass
#
#
# def list_part(request):
#     pass


# 编写返回数据的函数，以json的格式进行数据返回
# def info_list_vue(request):
#     if request.method == 'GET':
#         page = request.GET.get("page")
#         page_size = request.GET.get("page_size")
#         once_page = 10
#         if not page:
#             page = 1
#         if not page_size:
#             page_size = once_page
#         page = int(page)
#         page_size = int(page_size)
#         s_n = page / once_page
#         if page % once_page == 0:
#             s_n = int(page / once_page)
#             r_n = once_page
#         else:
#             s_n = int(page / once_page) + 1
#             r_n = once_page
#
#         total = pegRNAInfo.objects.all()  # 查询所有数据
#
#         select_start = (s_n - 1) * once_page * page_size
#         select_end = s_n * once_page * page_size
#         select_range = total[select_start:select_end]
#         return_start = (r_n - 1) * page_size
#         return_end = r_n * page_size
#         return_range = select_range[return_start:return_end]
#         if page <= 3:
#             page_range = [1, 2, 3, 4, 5]
#         else:
#             page_range = range(page - 2, page + 3)
#         return_list = []
#         for data in return_range:
#             return_list.append({
#                 "PegRNAID": data.PegRNAID,
#                 "Target": data.Target,
#                 "Spacer": data.Spacer,
#                 "PBS": data.PBS,
#                 "RT": data.RT,
#                 "NickingCoordinate": data.NickingCoordinate,
#
#             })
#         result = {"data": return_list, "page_range": ','.join([str(i) for i in page_range])}
#         return JsonResponse(result)

# return render(request, 'data_select_vue.html', {'result':result})
# def infoList_vue(request):
#     return render(request,'data_select_vue.html.html',locals())

# 文章查找
def getWhere(request):
    where = dict()
    article_title = request.POST.get('article_title', '')
    article_is_recommend = request.POST.get('article_is_recommend', '')
    member_id = request.POST.get('member_id', '')
    article_auditor = request.POST.get('article_auditor', '')
    if article_title:
        where['article_title__icontains'] = article_title
    if article_is_recommend != '':
        where['article_is_recommend'] = article_is_recommend
    if member_id:
        where['member'] = member_id
    if article_auditor != '':
        where['article_auditor'] = article_auditor
    print(where)
    return where

#
# def pegRNA_list(request):
#     result_list = pegRNAInfo.objects.order_by("PegRNAID")
#     num = len(result_list)
#     currentPage = int(request.GET.get('page', 1))
#     paginator = Paginator(result_list, 10)
#
#     if paginator.num_pages > 11:
#         if currentPage - 5 < 1:
#             pageRange = range(1, 11)
#         elif currentPage + 5 > paginator.num_pages:
#             pageRange = range(currentPage - 5, paginator.num_pages + 1)
#         else:
#             pageRange = range(currentPage - 5, currentPage + 5)
#     else:
#         pageRange = range(1, paginator.num_pages + 1)
#
#     try:
#         result_list = paginator.page(currentPage)
#     except PageNotAnInteger:
#         result_list = paginator.page(1)
#     except EmptyPage:
#         result_list = paginator.page(paginator.num_pages)
#
#     return render(request, 'data_select.html', locals())

# def info_list(request):
#     member_list = pegRNAInfo.objects.all()
#     return render(request,'info/info_list.html',locals())


# def list_part(request):
#     return render(request)

# def peg_list(request):
#     #
#     peg_list = pegRNAInfo.objects.order_by("-PegRNAID")
#     # 分页功能，一页7条数据
#     paginator = Paginator(peg_list, 7)
#     if request.method == 'GET':
#         # 默认显示第一页的数据
#         pegs = paginator.page(1)
#         return render(request, 'data_select.html', {'pegs': pegs})
#     # Ajax数据交互
#     if request.is_ajax():
#         page = request.GET.get('page')
#         try:
#             pegs = paginator.page(page)
#         # 如果页数不是整数，返回第一页
#         except PageNotAnInteger:
#             pegs = paginator.page(1)
#         # 如果页数不存在/不合法，返回最后一页
#         except InvalidPage:
#             pegs = paginator.page(paginator.num_pages)
#         peg_li = list(pegs.object_list.values())
#         # 分别为是否有上一页false/true，是否有下一页false/true，总共多少页，当前页面的数据
#         result = {'has_previous': pegs.has_previous(),
#                   'has_next': pegs.has_next(),
#                   'num_pages': pegs.paginator.num_pages,
#                   'peg_li': peg_li}
#         return render(request)
# import pandas as pd
# df=pd.
