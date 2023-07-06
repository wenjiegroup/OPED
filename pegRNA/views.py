import re

import pandas as pd
from django.shortcuts import render
from datetime import datetime
from pegRNA.utils import *
from pegRNA.models import *
from django.http import FileResponse
from django.core.paginator import Paginator, PageNotAnInteger, EmptyPage


def index(request):
    return render(request, "index.html")


def page_sequence(request):
    if request.method == "POST":
        input_dict = request.POST.dict()
        t = datetime.now().strftime('%y%m%d%H%M%S_')
        uid = get_uid(t)
        input_dict['UID'] = uid

        if len(request.FILES) == 0:
            info_list, err = run_by_seq(input_dict)
            if err == "NO ERROR":
                info_list['EditingScore'] = info_list.apply(decimal2, axis=1, args=('EditingScore',))

            start = 1
            return render(request, "Sequence.html", {"info_list": info_list, "err": err, "start": start, 'uid': uid})
        else:
            with open(f'Temp/User.{uid}.fa', 'wb+') as destination:
                for chunk in request.FILES['fa_in'].chunks():
                    destination.write(chunk)
            # form = UploadFileForm(request.POST, request.FILES)
            info_list, errs = run_by_seq_fa(input_dict)
            try:
                info_list['EditingScore'] = info_list.apply(decimal2, axis=1, args=('EditingScore',))
            except:
                pass
            start = 2
            return render(request, "Sequence.html",
                          {"info_list": info_list, "errs": errs, "start": start, 'uid': uid})
    return render(request, "Sequence.html")


def page_position(request):
    if request.method == "POST":
        input_dict = request.POST.dict()
        t = datetime.now().strftime('%y%m%d%H%M%S_')
        uid = get_uid(t)
        input_dict['UID'] = uid

        if len(request.FILES) == 0:
            # Whether or not an error occurs is encapsulated in the function of hjs_file_pos
            info_list, err = run_by_pos(input_dict)
            if err == "NO ERROR":
                info_list['EditingScore'] = info_list.apply(decimal2, axis=1, args=('EditingScore',))

            start = 1
            return render(request, "Position.html", {"info_list": info_list, 'err': err, 'start': start, 'uid': uid})
        else:
            with open(f'Temp/User.{uid}.bed', 'wb+') as destination:
                for chunk in request.FILES['bed_in'].chunks():
                    destination.write(chunk)
            # form = UploadFileForm(request.POST, request.FILES)
            info_list, errs = run_by_pos_bed(input_dict)
            try:
                info_list['EditingScore'] = info_list.apply(decimal2, axis=1, args=('EditingScore',))
            except:
                pass
            start = 2
            return render(request, "Position.html",
                          {"info_list": info_list, 'errs': errs, 'start': start, 'uid': uid})
    return render(request, "Position.html")


def page_opedvar(request):
    if request.method == "POST":
        input_dict = request.POST.dict()
        t = datetime.now().strftime('%y%m%d%H%M%S_')
        uid = get_uid(t)
        input_dict['UID'] = uid

        # print(input_dict)
        info_list, err = run_by_opedvar(input_dict)

        query = dict(queryType=input_dict['queryType'], queryItem=input_dict['queryItem'], PAM=input_dict['PAM'])

        # print(info_list)

        start = 1
        return render(request, "OPEDVar.html", {"info_list": info_list, "err": err, "start": start, 'uid': uid, "query": query})
    return render(request, "OPEDVar.html")


def history_page(request):
    p = request.path
    section = re.search(r'/([^/]+)/history', p).group(1)
    all_recs = eval(f"cookie_{section.lower()}_index").objects.all().order_by('-UID')

    # print(all_recs)
    if all_recs:
        paginator = Paginator(all_recs, 10)
        page = request.GET.get('page')
        try:
            contacts = paginator.page(page)
            # print(contacts)
        except PageNotAnInteger:
            page = 1
            contacts = paginator.page(page)
        except EmptyPage:
            page = paginator.num_pages
            contacts = paginator.page(page)

        page = int(page)

        pgrange = range(max(1, page - 4), min(paginator.num_pages + 1, page + 5))
        start_history = 1
        return render(request, f'{section}.html', {'contacts': contacts, 'start_history': start_history, 'pgrange': pgrange})
    # else:
    #     start_history = 1
    #     history_updateinfo = 'No History'
    #     return render(request, f'{section}.html', {"history": history_updateinfo})


def download_page(request):
    p = request.path
    section = re.search(r'/([^/]+)/download', p).group(1)
    uid = request.GET.get("uid", None)

    qrs = unique_id.objects.filter(UID__iregex=f"^{uid}")
    res0 = pd.DataFrame(None)
    for qr in qrs:
        res = eval(f"qr.cookie_{section.lower()}_set").all()
        res = pd.DataFrame(list(res.values()))
        res0 = pd.concat([res0, res])
    # filter_model = eval(f"rec.cookie_{section.lower()}_set").all()
    # filter_model = eval(f"cookie_{section.lower()}").objects.filter(UID=uid)
    # filter_model.df = pd.DataFrame(list(filter_model.values()))

    res0.to_csv(f"Temp/{section}_{uid}.csv")
    response = FileResponse(open(f"Temp/{section}_{uid}.csv", 'rb'))
    response['Content-Type'] = 'application/octet-stream'
    response['Content-Disposition'] = f'attachment;filename="{section}_{uid}.csv"'
    return response


def detail_page(request):
    p = request.path
    section = re.search(r'/([^/]+)/history/detail', p)
    section = section.group(1)
    uid = request.GET.get("uid", None)

    rec = unique_id.objects.get(UID=uid)
    filter_model = eval(f"rec.cookie_{section.lower()}_set").all()
    filter_model.df = pd.DataFrame(list(filter_model.values()))
    start_history = 2
    if section == 'OPEDVar':
        query = eval(f"rec.cookie_{section.lower()}_index")
        return render(request, f'{section}.html',
                      {'contacts_df': filter_model.df, 'start_history': start_history, 'uid': uid, 'query': query})

    # print(filter_model.df['Target3'])
    return render(request, f'{section}.html',
                  {'contacts_df': filter_model.df, 'start_history': start_history, 'uid': uid})


def page_help(request):
    pass
    return render(request, 'help.html')


def page_about(request):
    pass
    return render(request, 'about.html')

