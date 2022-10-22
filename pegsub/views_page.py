from django.shortcuts import render
from django.views.generic.base import View
from django.db.models import Q
from pure_pagination import Paginator, EmptyPage, PageNotAnInteger
import re
import os
import sys
import subprocess

class SearchView(View):
    def get(self, request):
        Chr_pos = request.GET.get("Chr", "")
        Specie = request.GET.get("Specie", "")
        Genename = request.GET.get("Genename", "")

        if Chr_pos != '':
            Chr = Chr_pos.split(':')[0]
            pos = Chr_pos.split(':')[1]
            start_pos = int(pos.split('-')[0])
            end_pos = int(pos.split('-')[1])

        if Specie == 'Human':
            if Chr_pos != '':
                RNAedits = Human_Rnaedit.objects.filter(Q(Chr=Chr),Q(position__range=[start_pos,end_pos]))
                if Genename:
                    RNAedits = RNAedits.filter(Genename__icontains=Genename)
            elif Genename != '':
                RNAedits = Human_Rnaedit.objects.filter(Genename__icontains=Genename)
            else:
                return render(request, "search.html", {
                "Error_message":"Please pass chr or gene arg",
                })
        else:
            return render(request, "search.html", {
            "Error_message":"Please pass specie arg",
            })

        AAChange_list = request.GET.get('AAChange',"")
        AAChange_str = AAChange_list
        if AAChange_list:
            AAChange_list = AAChange_list.split(',')
            RNAedits = RNAedits.filter(AAChange__in = AAChange_list)
            for key in AAChangeTypes:
                if key in AAChange_list:
                    AAChangeTypes[key] = 1
                else:
                    AAChangeTypes[key] = 0
        else:
            for key in AAChangeTypes:
                if key in AAChange_list:
                    AAChangeTypes[key] = 1
                else:
                    AAChangeTypes[key] = 0

        GeneRegion_list = request.GET.get('GeneRegion', "")
        GeneRegion_str = GeneRegion_list
        if GeneRegion_list:
            GeneRegion_list = GeneRegion_list.split(',')
            RNAedits = RNAedits.filter(GeneRegion__in = GeneRegion_list)
            for key in GeneRegionTypes:
                if key in GeneRegion_list:
                    GeneRegionTypes[key] = 1
                else:
                    GeneRegionTypes[key] = 0
        else:
            for key in GeneRegionTypes:
                if key in GeneRegion_list:
                    GeneRegionTypes[key] = 1
                else:
                    GeneRegionTypes[key] = 0


        Repeat_list = request.GET.get('Location', "")
        Repeat_str = Repeat_list
        if Repeat_list:
            Repeat_list = Repeat_list.split(',')
            RNAedits = RNAedits.filter(Repeat__in = Repeat_list)
            for key in RepeatTypes:
                if key in Repeat_list:
                    RepeatTypes[key] = 1
                else:
                    RepeatTypes[key] = 0
        else:
            for key in RepeatTypes:
                if key in Repeat_list:
                    RepeatTypes[key] = 1
                else:
                    RepeatTypes[key] = 0

        Method_list = request.GET.get('Method', "")
        Method_str = Method_list
        if Method_list:
            Method_list = Method_list.split(',')
            for key in MethodTypes:
                if key in Method_list:
                    MethodTypes[key] = 1
                else:
                    MethodTypes[key] = 0
            Method_list = [Methoddict[i] for i in Method_list]
            condition = Q()
            for i in Method_list:
                condition |= Q(method__icontains = i)
            RNAedits = RNAedits.filter(condition)

        else:
            for key in MethodTypes:
                if key in Method_list:
                    MethodTypes[key] = 1
                else:
                    MethodTypes[key] = 0

        Project_list = request.GET.get('Project', "")
        Project_str = Project_list
        if Project_list:
            Project_list = ';'.join(Project_list.split(','))
            RNAedits = RNAedits.filter(Project__icontains = Project_list)
            for key in Human_Project:
                if key in Project_list:
                    Human_Project[key] = 1
                else:
                    Human_Project[key] = 0
        else:
            for key in Human_Project:
                if key in Project_list:
                    Human_Project[key] = 1
                else:
                    Human_Project[key] = 0

        dbSNP_list = request.GET.get('dbSNP', "")
        dbSNP_str = dbSNP_list
        if dbSNP_list:
            dbSNP_list = dbSNP_list.split(',')
            if 'yes' in dbSNP_list and 'no' not in dbSNP_list:
                RNAedits = RNAedits.exclude(dbSNP = '-')
            elif 'yes' not in dbSNP_list and 'no' in dbSNP_list:
                RNAedits = RNAedits.filter(dbSNP = '-')
            for key in dbSNPTypes:
                if key in dbSNP_list:
                    dbSNPTypes[key] = 1
                else:
                    dbSNPTypes[key] = 0
        else:
            for key in dbSNPTypes:
                if key in dbSNP_list:
                    dbSNPTypes[key] = 1
                else:
                    dbSNPTypes[key] = 0

        Number = RNAedits.count()
        try:
            page = request.GET.get('page', 1)
        except PageNotAnInteger:
            page = 1

        p = Paginator(RNAedits, 10, request=request)

        RNAedits = p.page(page)
        return render(request, "result.html", {
            "Chr":Chr_pos,
            "Specie":Specie,
            "Genename":Genename,
            "RNAedits":RNAedits,
            "all_AAChanges":AAChangeTypes,
            "all_locations":RepeatTypes,
            "all_generegions":GeneRegionTypes,
            "all_projects":Human_Project,
            "all_methods":MethodTypes,
            "all_dbsnp":dbSNPTypes,
            "number":Number,
            "AAChange_str":AAChange_str,
            "GeneRegion_str":GeneRegion_str,
            "Repeat_str":Repeat_str,
            "Method_str":Method_str,
            "Project_str":Project_str,
            "dbSNP_str":dbSNP_str,
            })
