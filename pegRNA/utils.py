import random
import subprocess
import torch
import os
import re
from django.db.utils import IntegrityError as ie
import pandas as pd
from decimal import Decimal
from Bio import SeqIO
from pegRNA_PredictingCodes import evaluate_model, train_model
from pegRNA_PredictingCodes.predict_efficiency_of_ClinVar import read_data_of_ClinVar_file
from pegRNA.models import *


# get .xx  Decimal
def decimal2(x, i):  # i = 'score'
    return Decimal(x[i]).quantize(Decimal("0.00"))


def random_str(n=4):
    s = ''
    base_str = 'abcdefghijklmnopqrstuvwxyz0123456789'
    l = len(base_str) - 1
    for i in range(n):
        s += base_str[random.randint(0, l)]
    return s


def get_uid(t):
    try:
        s = t + random_str()  # 230601120156_ax1d
        uid = unique_id.objects.create(UID=s)
    except ie:
        uid = get_uid(t)
    return uid


# save dict / dicts into models (database)
def rec2model(bulk, rec, model, **kwargs):
    fileds = [i.name for i in model._meta.get_fields()]
    fileds.remove('id')
    if bulk:  # True
        temp = []
        for i in rec:
            x = rec[i]
            for j in kwargs:
                x[j] = kwargs[j]
            x2 = {j: x[j] for j in fileds}
            temp.append(model(**x2))
        model.objects.bulk_create(temp)
    else:  # False
        x = rec
        for j in kwargs:
            x[j] = kwargs[j]
        x2 = {j: x[j] for j in fileds}
        model.objects.create(**x2)

    return


def run_model(input_dict, tp):
    # scripts
    scps = dict(Sequence="scripts/A0_PE_library_for_user_sequence.pl",
                Position="scripts/A1_PE_library_for_user_position.pl",
                OPEDVar="scripts/A2_search_database.pl")
    # parameters to be used in scripts
    keys = dict(Sequence=['SEQUENCE', 'PAM', 'CUT_SIZE', 'MIN_DISGRNA', 'MAX_DISGRNA', 'MAX_PBS', 'MAX_RT',
                          'MIN_PBS', 'MIN_RT', 'HOMOLOGY'],
                Position=['Chromosome', 'Position', 'Pattern', 'PAM', 'CUT_SIZE', 'MIN_DISGRNA', 'MAX_DISGRNA',
                          'MAX_PBS', 'MAX_RT', 'MIN_PBS', 'MIN_RT', 'HOMOLOGY'],
                OPEDVar=['queryType', 'queryItem', 'PAM', 'DIRECTION'])

    sub = subprocess.Popen(
        ["perl", scps[tp]] + [input_dict[i] for i in keys[tp]] + [input_dict['UID'].UID],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf8'
    )
    sub.wait()
    outerror = sub.stderr.read()
    return outerror


def run_by_seq(seq_dict):
    seq_dict['SEQUENCE'] = ''.join(re.split(r'\s+', seq_dict['SEQUENCE']))
    # seq_dict['SEQUENCE'].strip()
    if seq_dict['SEQUENCE'] == '':
        err = "Please enter a sequence!"  # xuhao
        o = dict(Target=err, PBS=err, RT=err, max=err)
        return o, err

    # whether error appear,we catch the error for
    err = run_model(seq_dict, 'Sequence')
    if err == "":
        err = "NO ERROR"

        device = torch.device(f"cpu")
        transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                             model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')  # load model
        df = pd.read_table(f'Temp/Sequence.request.User.{seq_dict["UID"]}.txt')
        os.remove(f'Temp/Sequence.request.User.{seq_dict["UID"]}.txt')
        df = df.dropna()

        if df.shape[0] == 0:
            err = "Sorry... there is no proper candidates found, please check your sequence"  # xuhao
            o = dict(Target=err, PBS=err, RT=err, max=err)
        else:
            # when the shape of the df is not equal to zero means that the shell run in normal
            # so the error
            data_X = read_data_of_ClinVar_file(df)
            o = evaluate_model.transformer_predictor_order3_file_update(transformer, data_X, 1024, device,
                                                                        seq_dict['TOP_N'], df)
            o = o.replace('Na', None)
            o2 = o.to_dict(orient='index')
            rec2model(True, o2, cookie_sequence, FID=seq_dict['UID'])
            rec2model(False, seq_dict, cookie_sequence_index)
    else:
        o = dict(Target=err, PBS=err, RT=err, max=err)
    return o, err


def run_by_seq_fa(seq_dict):
    seqs = [seq for seq in SeqIO.parse(f"Temp/User.{seq_dict['UID']}.fa", "fasta")]
    outputs = pd.DataFrame()
    errs = []
    uid0 = seq_dict['UID'].UID
    for seq in seqs:
        seq_dict['UID'] = get_uid(uid0)  # new UID
        seq_dict['SEQUENCE'] = str(seq.seq)

        if seq_dict['SEQUENCE'] == '':
            err = "Please enter a sequence!"  # xuhao
        else:
            err = run_model(seq_dict, 'Sequence')
            if err == "":
                device = torch.device(f"cpu")
                # load model
                transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                                     model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')
                df = pd.read_table(f'Temp/Sequence.request.User.{seq_dict["UID"]}.txt')
                os.remove(f'Temp/Sequence.request.User.{seq_dict["UID"]}.txt')
                df = df.dropna()

                if df.shape[0] == 0:
                    err = "Sorry... there is no proper candidates found, please check your sequence"  # xuhao
                else:
                    err = "NO ERROR"
                    data_X = read_data_of_ClinVar_file(df)
                    o = evaluate_model.transformer_predictor_order3_file_update(transformer, data_X, 1024, device,
                                                                                seq_dict['TOP_N'], df)
                    o = o.replace('Na', None)
                    o['SequenceName'] = seq.name
                    o2 = o.to_dict(orient='index')
                    rec2model(True, o2, cookie_sequence, FID=seq_dict['UID'])
                    rec2model(False, seq_dict, cookie_sequence_index)

                    outputs = pd.concat([outputs, o])

        errs.append([seq.name, err])
    return outputs, errs


def run_by_pos(pos_dict):
    err = run_model(pos_dict, 'Position')
    if err == "":
        err = "NO ERROR"

        device = torch.device(f"cpu")
        transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                             model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')  # load model
        df = pd.read_table(f'Temp/Position.request.User.{pos_dict["UID"]}.txt', header=0)
        os.remove(f'Temp/Position.request.User.{pos_dict["UID"]}.txt')
        df = df.dropna()
        if df.shape[0] == 0:
            err = "Sorry... there is no proper candidates found, please check your sequence"  # xuhao
            o = {'Strand': None, 'Spacer': None, 'PAM ': None,
                 'PBS': None, 'RT': None, 'EditToNickDistance': None,
                 'sgRNAPAM': None, 'sgRNASpacer': None, 'EditingScore': None, 'time': None
                 }
        else:
            data_X = read_data_of_ClinVar_file(df)
            o = evaluate_model.transformer_predictor_order3_file_pos(transformer, data_X, 1024, device,
                                                                     pos_dict['TOP_N'], df)
            o = o.replace('Na', None)
            o2 = o.to_dict(orient='index')
            rec2model(True, o2, cookie_position, FID=pos_dict['UID'])
            rec2model(False, pos_dict, cookie_position_index)
    else:
        o = dict(Target=err, PBS=err, RT=err, max=err)
    return o, err


def run_by_pos_bed(pos_dict):
    poses = []
    with open(f'Temp/User.{pos_dict["UID"]}.bed', 'r') as bed:
        r = bed.readlines()
        for l in r:
            if l.startswith('#'):
                continue
            pos = l.strip().split('\t')
            if len(pos) == 5:
                poses.append(pos)
    outputs = pd.DataFrame()
    errs = []
    uid0 = pos_dict['UID'].UID
    for pos in poses:
        pos_dict['UID'] = get_uid(uid0)  # new UID
        pos_name = pos[4]
        pos_dict['Chromosome'] = pos[0]  # .lower()
        pos_dict['Position'] = pos[2]
        pos_dict['Pattern'] = pos[3]

        err = run_model(pos_dict, 'Position')
        if err == "":
            device = torch.device(f"cpu")
            transformer = train_model.load_model(device, model_dir='pegRNA_PredictingCodes/Model_Trained',
                                                 model_name='pegRNA_Model_Merged_saved.order3_decoder.pt')  # load model
            df = pd.read_table(f'Temp/Position.request.User.{pos_dict["UID"]}.txt', header=0)
            os.remove(f'Temp/Position.request.User.{pos_dict["UID"]}.txt')
            df = df.dropna()
            if df.shape[0] == 0:
                err = "Sorry... there is no proper candidates found, please check your sequence"  # xuhao
            else:
                err = "NO ERROR"
                data_X = read_data_of_ClinVar_file(df)
                o = evaluate_model.transformer_predictor_order3_file_pos(transformer, data_X, 1024, device,
                                                                         pos_dict['TOP_N'], df)
                o = o.replace('Na', None)
                o['PositionName'] = pos_name
                o2 = o.to_dict(orient='index')
                rec2model(True, o2, cookie_position, FID=pos_dict['UID'])
                rec2model(False, pos_dict, cookie_position_index)

                outputs = pd.concat([outputs, o])

        errs.append([pos_name, err])
    return outputs, errs


def run_by_opedvar(opedvar_dict):
    err = run_model(opedvar_dict, 'OPEDVar')
    if err:
        df = pd.DataFrame(
            columns=['AlleleID', 'Type', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'Spacer',
                     'PBS', 'RT', 'PAM'])
    else:
        err = "NO ERROR"
        df = pd.read_table(f'Temp/OPEDVar.request.{opedvar_dict["UID"]}.txt', header=0)
        os.remove(f'Temp/OPEDVar.request.{opedvar_dict["UID"]}.txt')
        if df.shape[0] == 0:
            err = "Sorry... there is no proper candidates found, please try anothor one"  # xuhao
            df = {'Strand': None, 'Spacer': None, 'PAM ': None,
                 'PBS': None, 'RT': None, 'EditToNickDistance': None,
                 'sgRNAPAM': None, 'sgRNASpacer': None, 'EditingScore': None, 'time': None
                 }
        else:
            df = df.replace('Na', None)

            df2 = df.to_dict(orient='index')
            rec2model(True, df2, cookie_opedvar, FID=opedvar_dict['UID'])
            rec2model(False, opedvar_dict, cookie_opedvar_index)

    return df, err



