from pip._internal.utils import logging

from pegRNA_PredictingCodes import evaluate_model
import torch
import pandas as pd
from pegRNA_PredictingCodes import train_model
from pegRNA_PredictingCodes import predict_efficiency_of_ClinVar
import sys
import subprocess
from threading import main_thread

from pegRNA_PredictingCodes.predict_efficiency_of_ClinVar import read_data_of_ClinVar, read_data_of_ClinVar_file, \
    read_data_of_Single
import time
from pegsub.models import *

sys.path.append("D:/Pycharm_workplace/pegRNA")
import sys
import os


# seq = "CGCCGCCGCCTCTACTGGGGCTTCTTCTCGGGCCCGGCCGCGTCAAGCCGGGGGGGCGCTGGCCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC"
# seq = "CGCCGCCGCCTGGGGCTTCTTCTCGGGCCCGGCCGCGTCAAGCCGGGGGGGCGCTGGCCGA(/T)GGCCGCCTGGCAACTCTGCGACTACTACCTGCC"
'''hjs_file暂时不使用
def hjs_file(seq_input):
    hjs_seq_to_file(seq_input)

    device = torch.device(f"cuda:{1}" if torch.cuda.is_available() else "cpu")
    transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                         model_name='pegRNA_Model_Merged_saved.order3_decoder_ori.pt')  # load model

    # df = pd.read_table(f'Sequence/pegRNA.User.txt', header=0)
    #新的update
    df=pd.read_table(f'Sequence/pegRNA_April.User.txt', header=0)

    ####
    df = df.dropna()
    print(df)
    data_X = read_data_of_ClinVar_file(df)

    o = evaluate_model.transformer_predictor_order3_file(transformer, data_X, 1024, device)

    print(o)

    # seq_input=
    return o

'''
# 颜色展示
def get_Target1():
    pass


##data_hjs_update
def hjs_file_update(seq_dict, time_str):
    # whether error appear,we catch the error for
    ERROR = hjs_seq_to_file_update(seq_dict)
    print("meiyou")
    if ERROR=="":
        ERROR="None ERROR"
    print(ERROR)
    print("meiyouend")

    device = torch.device(f"cpu")
    transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                         model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')  # load model
    # df = pd.read_table(f'Sequence/pegRNA_update.User.txt', header=0)
    df=pd.read_table(f'Sequence/pegRNA_April_update.User.txt')
    df = df.dropna()

    if df.shape[0] == 0 :
        # when the shape of the df is zero means that the shell of perl occur error.
        # to avoid the webserver interrupt we catch the error and render to the html showing for the error
        # the error invalidation was encapsulater via perl shell from @liufeng

        o = {
            'Target': ERROR,
            'PBS': ERROR,
            'RT': ERROR,
            'max': ERROR
        }

        return o, ERROR
    else:
        # when the shape of the df is not equal to zero means that the shell run in normal
        # so the error
        print(df)
        data_X = read_data_of_ClinVar_file(df)
        o = evaluate_model.transformer_predictor_order3_file_update(transformer, data_X, 1024, device,
                                                                    seq_dict['top_n'])
        # o.to_csv(f"Sequence/show_pegRNA_update.User.csv")
        o.to_csv(f"Sequence/show_pegRNA_update_April.csv")

        '''此部分注释信息为四月份更新前 sequence 即update部分的输出内容cookie信息
        for index, row in o.iterrows():
            temp_cookie = cookie_sequence()
            temp_cookie.Target = row['Target']
            temp_cookie.PBS = row['PBS']
            temp_cookie.RT = row['RT']
            temp_cookie.Efficiency = row['score']
            temp_cookie.time = time_str
            temp_cookie.save()
        '''
        #此部分为更新后update部分的输出内容的cookie信息
        for index,row in o.iterrows():
            temp_cookie=cookie_sequence_April()
            temp_cookie.Strand=row['Strand']
            temp_cookie.Spacer=row['Spacer']
            temp_cookie.PAM=row['PAM']
            temp_cookie.PBS=row['PBS']
            temp_cookie.RT=row['RT']
            temp_cookie.EditToNickDistance=row['EditToNickDistance']
            temp_cookie.sgRNASpacer=row['sgRNASpacer']
            temp_cookie.sgRNAPAM=row['sgRNAPAM']
            temp_cookie.NickToNickDistance=row['NickToNickDistance']
            temp_cookie.Efficiency=row['score']
            temp_cookie.time=time_str
            temp_cookie.save()


        '''此部分为update部分更新前的
        temp_cookie_index = cookie_sequence_index()
        temp_cookie_index.PAM = seq_dict['pam_type']
        temp_cookie_index.CUT_SIZE = seq_dict['dis_nickase']
        temp_cookie_index.MAX_PBS = seq_dict['max_len_PBS']
        temp_cookie_index.MAX_RT = seq_dict['max_len_RT']
        temp_cookie_index.MIN_PBS = seq_dict['min_len_PBS']
        temp_cookie_index.MIN_RT = seq_dict['min_len_RT']
        temp_cookie_index.number_show = seq_dict['top_n']
        temp_cookie_index.input_sequence = seq_dict['user_seq']
        temp_cookie_index.time = time_str
        temp_cookie_index.save()
'''

        #此部分为更新后的
        temp_cookie_index=cookie_sequence_April_index()
        temp_cookie_index.PAM=seq_dict['pam_type']
        temp_cookie_index.CUT_SIZE=seq_dict['dis_nickase']
        temp_cookie_index.MAX_PBS = seq_dict['max_len_PBS']
        temp_cookie_index.MAX_RT = seq_dict['max_len_RT']
        temp_cookie_index.MIN_PBS = seq_dict['min_len_PBS']
        temp_cookie_index.MIN_RT = seq_dict['min_len_RT']
        temp_cookie_index.number_show = seq_dict['top_n']

        temp_cookie_index.MIN_DISGRNA=seq_dict['min_dis_sgRNA_nickase']
        temp_cookie_index.MAX_DISGRNA=seq_dict['max_dis_sgRNA_nickase']

        temp_cookie_index.HOMOLOGY=seq_dict['homology']
        temp_cookie_index.input_sequence=seq_dict['user_seq']
        temp_cookie_index.time=time_str
        temp_cookie_index.save()



        return o, ERROR


##data_hjs_pos
def hjs_file_pos(pos_input, time_str):
    ERROR = hjs_seq_to_file_pos(pos_input)
    if ERROR=="":
        ERROR="None ERROR"

    device = torch.device(f"cpu")
    transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                         model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')  # load model
    try:
        # 四月更新前的df存储地点
        # df = pd.read_table(f'Sequence/pegRNA_pos.User.txt', header=0)
        #四月更新后
        df=pd.read_table(f'Position/pegRNA_April_pos.User.txt', header=0)
        print(df)
    except IOError:
        print("ERROR:MEIYOUZHAODAO")
    df = df.dropna()
    print(df.shape[0])
    if (df.shape[0] == 0):
        print(df.shape[0])
        ERROR="THIS POSITION DIDN'T HAVE ACTG"
        o= { 'Strand':None, 'Spacer':None, 'PAM ':None,
            'PBS':None, 'RT' : None, 'EditToNickDistance' :None,
        'sgRNAPAM' : None, 'sgRNASpacer' : None, 'Efficiency' : None, 'time' : None
        }
    else:
        data_X = read_data_of_ClinVar_file(df)
        o = evaluate_model.transformer_predictor_order3_file_pos(transformer, data_X, 1024, device, pos_input['top_n'])
        # o.to_csv(f"Sequence/show_pegRNA_position.User.csv") 四月更改前
        o.to_csv(f"Position/show_pegRNA_position_April.User.csv")

        for index, row in o.iterrows():
            '''
            cookie_position.objects.create(Target=row['Target'], PBS=row['PBS'], RT=row['RT'],
                                           #若要存入数据库的话用Efficiency=Decimal(row['score']).quantize(Decimal("0.00"))
                                           Efficiency=row['score'], time=time_str)
            '''
            cookie_position_April.objects.create(Strand=row['Strand'],Spacer=row['Spacer'],PAM=row['PAM'],
                                                 PBS=row['PBS'],RT=row['RT'],EditToNickDistance=row['EditToNickDistance'],
                                                 sgRNAPAM=row['sgRNAPAM'],sgRNASpacer=row['sgRNASpacer'],NickToNickDistance=row['NickToNickDistance'],Efficiency=row['score'],time=time_str)
        '''此部分为更新前cookie更新
        cookie_position_index.objects.create(PAM=pos_input['pam_type'], CUT_SIZE=pos_input['dis_nickase'],
                                             MAX_PBS=pos_input['max_len_PBS'], MAX_RT=pos_input['max_len_RT'],
                                             MIN_PBS=pos_input['min_len_PBS'], MIN_RT=pos_input['min_len_RT'],
                                             ###
                                             #数据库信息得更新
                                             number_show=pos_input['top_n'], input_chr=pos_input['chr'],
                                             input_position=pos_input['start_position'],
                                             input_pattern=pos_input['edit_pattern'],
                                             time=time_str
                                             )
        '''
        cookie_position_April_index.objects.create(PAM=pos_input['pam_type'],CUT_SIZE=pos_input['dis_nickase'],
                                                   MAX_PBS=pos_input['max_len_PBS'],MAX_RT=pos_input['max_len_RT'],
                                                   MIN_PBS=pos_input['min_len_PBS'],MIN_RT=pos_input['min_len_RT'],
                                                   number_show=pos_input['top_n'],

                                                   MIN_DISGRNA=pos_input['min_dis_sgRNA_nickase'],
                                                   MAX_DISGRNA=pos_input['max_dis_sgRNA_nickase'],

                                                   HOMOLOGY=pos_input['homology'],

                                                   input_chr=pos_input['chr'],
                                                   input_position=pos_input['start_position'],
                                                   input_pattern=pos_input['edit_pattern'],
                                                   time=time_str
                                                   )



    return o, ERROR


# test

def hjs_test(Target, PBS, RT):
    device = torch.device(f"cpu")
    data = {'Target(47bp)': [Target],
            'PBS': [PBS],
            'RT': [RT], }

    data = pd.DataFrame(data)
    # print(data_x)

    data_X = predict_efficiency_of_ClinVar.read_data_of_Single(data)

    data_X.to_csv(f'Data/temp.txt', header=1)
    print("hjs_change")
    transformer = train_model.load_model(device,
                                         model_dir='pegRNA_PredictingCodes/Model_Trained',
                                         model_name='pegRNA_Model_Merged_saved.order3_decoder_ori.pt', )

    o1 = evaluate_model.transformer_predictor_order3(transformer, data_X, 1024, device)
    print(o1)
    o = str(o1[0])

    return o


"""
1.根据输入序列生成文件
  ｛
    输入序列
    调用perl脚本
    生成一个文件  
  ｝
2.根据生成文件预测效率最高值并组合
  ｛
  方案一：
   1.创建一个字典，存储效率最高的TARGET\PBS\RT组合
   2.由于字典存储的组合都是数字，所以需要一个函数碱基互换的函数
   3.输出这个result
   方案二
   1.获取最大效率
   2.获取最大效率的 index
   3.返回文件index获取PBS RT tartget
  ｝
"""


# def predict_ClinVar_order3_file(gpu_id=1):
#     device = torch.device(f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu")
#     transformer = train_model.load_model(device, model_dir='Model_Trained',
#                                          model_name='pegRNA_Model_Merged_saved.order3_decoder_ori.pt')  # load model
#     data_X = read_data_of_ClinVar(type)
#
#     o = evaluate_model.transformer_predictor_order3(transformer, data_X, 1024, device)[0]
#
#     pd.DataFrame({'Efficiency': o}).to_csv(f'Output/Efficiency.pegRNA.GRCh38.{type}.txt',
#                                            sep='\t', index=False)


# def hjs_seq_to_file(seq_input):
#     perl_input = seq_input
#
#     subp = subprocess.Popen(["perl",
#                              # "A0_PE_library_for_user_sequence_update.pl",
#                              'A0_PE_library_for_user_sequence_April.pl',
#                              perl_input], stdin=subprocess.PIPE,
#                             stdout=subprocess.PIPE, encoding='utf8')
#     subp.wait()
#     print("Seq2")

    # subprocess.Popen.wait()


def hjs_seq_to_file_update(dict_seq):
    input = dict_seq
    # valation

    sub = subprocess.Popen(
        # ["perl", "A0_PE_library_for_user_sequence_update.pl",
        #  input['user_seq'],
        #  input['pam_type'],
        #  input['dis_nickase'],
        #  input['max_len_PBS'],
        #  input['max_len_RT'],
        #  input['min_len_PBS'],
        #  input['min_len_RT']
        #  ],
        ["perl","A0_PE_library_for_user_sequence_April.pl",
         input['user_seq'],
         input['pam_type'],
         input['dis_nickase'],
         #4.24 更新模型加字典
         input['min_dis_sgRNA_nickase'],
         input['max_dis_sgRNA_nickase'],
         ######
         input['max_len_PBS'],
         input['max_len_RT'],
         input['min_len_PBS'],
         input['min_len_RT'],
         # 更新模型 加字典 4.24
         input['homology'],
         ##
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf8'
    )
    sub.wait()
    outerror = sub.stderr.read()
    return outerror


# pos
def hjs_seq_to_file_pos(pos_seq):
    input = pos_seq

    sub = subprocess.Popen(
        ["perl", "A1_PE_library_for_user_position_April.pl",
         input['chr'],
         input['start_position'],
         input['edit_pattern'],
         input['pam_type'],
         input['dis_nickase'],
         #四月更新字段
         input['min_dis_sgRNA_nickase'],
         input['max_dis_sgRNA_nickase'],
         #四月更新字段end
         input['max_len_PBS'],
         input['max_len_RT'],
         input['min_len_PBS'],
         input['min_len_RT'],
         #更新模型
         input['homology']
         #end

         ],

        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8'

    )
    sub.wait()
    outerror = sub.stderr.read()
    return outerror


# search
def hjs_search(search_dict, time_str):
    input = search_dict
    # 删除上一次的结果
    # os.remove(f'User_Sequence/Database/output_of_search.txt')
    sub = subprocess.Popen(
        ["perl",
         # "A2_search_database.pl",
         "A2_search_database_April.pl",
         # input['alleleID'],
         input['queryType'],
         input['queryItem'],
         input['pam_type'],
         input['direction']
         ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf8'
    )
    sub.wait()

    outerror = sub.stderr.read()
    return outerror


def hjs_file_database(search_dict, time_str):
    ERROR = hjs_search(search_dict, time_str)
    if ERROR:
        df = pd.DataFrame(
            columns=['AlleleID', 'Type', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'Spacer',
                     'PBS', 'RT', 'PAM'])
        print(ERROR)
    else:
        df = pd.read_table(f'Data/output_of_search.txt', header=0)


        df['PAM']=search_dict['pam_type']
        print(df)
        # df.to_csv(f'Sequence/show_pegRNA_database.User.csv')
        df.to_csv(f'Data/show_pegRNA_database_April.User.csv')
        for index, row in df.iterrows():

            # cookie_databse.objects.create(AlleleID=row['AlleleID'], Type=row['Type'],
            #                               Chromosome=row['Chromosome'], Start=row['Start'],
            #                               Stop=row['Stop'], ReferenceAllele=row['ReferenceAllele'],
            #                               AlternateAllele=row['AlternateAllele'], Spacer=row['Spacer'],
            #                               PBS=row['PBS'], RT=row['RT'],
            #                               PAM=row['PAM'], time=time_str
            #                               )
            '''四月更新'''
            cookie_database_April.objects.create(AlleleID=row['AlleleID'], Type=row['Type'],
                                          GeneID=row['GeneID'],GeneSymbol=row['GeneSymbol'],
                                          HGNC_ID=row['HGNC_ID'],
                                          Chromosome=row['Chromosome'], Start=row['Start'],
                                          Stop=row['Stop'], ReferenceAllele=row['ReferenceAllele'],
                                          AlternateAllele=row['AlternateAllele'], Spacer=row['Spacer'],
                                          PBS=row['PBS'], RT=row['RT'],
                                          sgRNASpacer=row['sgRNASpacer'],
                                          EditToNickDistance=row['EditToickDistance'],
                                          NickToNickDistance=row['NickToNickDistance'],
                                          PAM=row['PAM'],
                                          time=time_str

            )
        # cookie_database_index.objects.create(ALLELEID=search_dict['alleleID'], PAM=search_dict['pam_type'],
        #                                          DIRECTION=search_dict['direction'], time=time_str
        #                                          )
        '''四月更新'''
        cookie_database_April_index.objects.create(queryType=search_dict['queryType'],queryItem=search_dict['queryItem'],
                                                   PAM=search_dict['pam_type'],DIRECTION=search_dict['direction'],
                                                   time=time_str
                                                   )
    return df,ERROR
    # if sub.stderr.read() is None:
    #     ERROR = "1"
    #     return df, ERROR
    # else:
    #     ERROR=sub.stderr.read()
    #     return df, ERROR
