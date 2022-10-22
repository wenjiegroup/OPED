from django.shortcuts import render
from pegsub import models


def request_update_index(request):
    # 从 URL 中取参数
    page_num = request.GET.get("page")
    print(page_num, type(page_num))  # page_num 为 str 类型
    # 书籍总数
    total_count = models.cookie_sequence_index.objects.all().order_by("-time").count()
    # 每一页显示多少条数据
    per_page = 10
    # 总共需要多少页码来显示
    total_page, m = divmod(total_count, per_page)
    # 如果还有数据
    if m:
        total_page += 1
    try:
        page_num = int(page_num)
        # 如果输入的页码数超过了最大的页码数，默认返回最后一页
        if page_num > total_page:
            page_num = total_page
        # 如果输入的页码数小于 1，则返回第一页
        if page_num < 1:
            page_num = 1
    except Exception as e:
        # 当输入的页码不是正经数字的时候 默认返回第一页的数据
        page_num = 1
    # 定义两个变量保存数据从哪儿取到哪儿
    data_start = (page_num - 1) * 10
    data_end = page_num * 10
    # 页面上最多展示的页码
    max_page = 11
    # 如果总页码数小于页面上最多展示的页码
    if total_page < max_page:
        max_page = total_page
    half_max_page = max_page // 2
    # 页面上展示的页码的开始页
    page_start = page_num - half_max_page
    # 页面上展示的页码的结束页
    page_end = page_num + half_max_page
    # 如果当前页减一半比 1 小
    if page_start <= 1:
        page_start = 1
        page_end = max_page
    # 如果当前页加一半比总页码还大
    if page_end > total_page:
        page_end = total_page
        page_start = total_page - max_page + 1
    all_book = models.cookie_sequence_index.objects.all().order_by("-time")[data_start:data_end]
    # 拼接 html 的分页代码
    html_list = []
    # 添加首页按钮
    html_list.append(
        '<li><a href="/sequence_index_list/?page=1" rel="external nofollow" rel="external nofollow" rel="external nofollow" >首页</a></li>')
    # 如果是第一页，就没有上一页
    if page_num <= 1:
        html_list.append(
            '<li class="disabled"><a href="#" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" ><span aria-hidden="true">«</span></a></li>'.format(
                page_num - 1))
    else:
        # 加一个上一页的标签
        html_list.append(
            '<li><a href="/sequence_index_list/?page={}" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" ><span aria-hidden="true">«</span></a></li>'.format(
                page_num - 1))
    # 展示的页码
    for i in range(page_start, page_end + 1):
        # 给当前页添加 active
        if i == page_num:
            tmp = '<li class="active"><a href="/sequence_index_list/?page={0}" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" >{0}</a></li>'.format(
                i)
        else:
            tmp = '<li><a href="/sequence_index_list/?page={0}" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" >{0}</a></li>'.format(
                i)
        html_list.append(tmp)
        # 如果是最后一页，就没有下一页
    if page_num >= total_page:
        html_list.append(
            '<li class="disabled"><a href="#" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" ><span aria-hidden="true">»</span></a></li>')
    else:
        html_list.append(
            '<li><a href="/sequence_index_list/?page={}" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" ><span aria-hidden="true">»</span></a></li>'.format(
                page_num + 1))
        # 添加尾页按钮
    html_list.append(
        '<li><a href="/sequence_index_list/?page={}" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" rel="external nofollow" >尾页</a></li>'.format(
            total_page))
    page_html = "".join(html_list)  # 拼接 html 的分页代码
    return render(request, "predict_file_update.html", {"books": all_book, "page_html": page_html})
