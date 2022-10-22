from django.shortcuts import render

from django.http import HttpResponse, FileResponse, StreamingHttpResponse, JsonResponse

from pegRNA_PredictingCodes.predict_efficiency_of_ClinVar import read_data_of_Single
from pegsub import models, hjs_change


def Simple_Version(request):
    if request.method =="method":
        return render(request,"Simple_version.html")
