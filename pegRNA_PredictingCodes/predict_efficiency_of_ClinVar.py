#!/usr/bin/env python3
import numpy as np
import pandas as pd
import torch

from pegRNA_PredictingCodes import train_model, evaluate_model
from django.db import models


def complement_seq(seq):
    """get the complementary sequence of the input sequence
    Args:
        seq: the input sequence
    Returns:
        complementary sequence
    """

    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([complement_map[s.upper()] for s in list(seq)])


def reverse_seq(seq):
    """get the reversed sequence of the input sequence
    Args:
        seq: the input sequence
    Returns:
        reversed sequence
    """

    return "".join(reversed(list(seq)))


def read_data_of_ClinVar(type='SNV'):
    """obtain the data in Article_293T for various types of regressions"""

    df = pd.read_table(f'Data/pegRNA.GRCh38.{type}.txt', header=0)
    # df = df.loc[:, ['spacer sequence', 'RT', 'PBS', 'Original (Bold: sgRNA)', 'Efficiency(%)']]
    # df = df.iloc[:10, :]
    data = {'Target': [], 'Target_o2': [], 'Target_o3': [], 'RT': [], 'RT_o2': [], 'RT_o3': [],
            'PBS': [], 'PBS_o2': [], 'PBS_o3': []}

    len_RT = []
    len_PBS = []

    for i in range(len(df)):
        temp = df.iloc[i]
        len_RT.append(len(temp['RT']))
        len_PBS.append(len(temp['PBS']))
    MAX_RT = max(len_RT)  # 31
    MAX_PBS = max(len_PBS)  # 17
    max_len_Target = 47
    print(f'Maximum length of (Target seq, PBS, RT): ({max_len_Target}, {MAX_PBS}, {MAX_RT})')

    id2char = list('ACGT')
    char2id = {char: i + 1 for i, char in enumerate(id2char)}
    char2id_o2 = {f'{char}{char_j}': i * len(id2char) + j + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char)}
    char2id_o3 = {f'{char}{char_j}{char_k}': i * len(id2char) * len(id2char) + j * len(id2char) + k + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char) for k, char_k in
                  enumerate(id2char)}
    for i, row in df.iterrows():
        seq = reverse_seq(complement_seq(row['PBS'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_PBS):
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_PBS - 1):
            temp.append(0)
        data['PBS_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_PBS - 2):
            temp.append(0)
        data['PBS_o3'].append(temp)

        seq = reverse_seq(complement_seq(row['RT'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_RT):
            temp.append(0)
        data['RT'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_RT - 1):
            temp.append(0)
        data['RT_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_RT - 2):
            temp.append(0)
        data['RT_o3'].append(temp)

        # seq = row['Original (Bold: sgRNA)'].upper()
        seq = row['Target(47bp)'].upper()
        # data['Target'].append([char2id[seq[j]] for j in range(0, len(seq))])
        # data['Target_o2'].append([char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)])
        # data['Target_o3'].append([char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)])
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_Target):
            temp.append(0)
        data['Target'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_Target - 1):
            temp.append(0)
        data['Target_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_Target - 2):
            temp.append(0)
        data['Target_o3'].append(temp)

        # data['Efficiency'].append(row['Efficiency(%)'] / 1)

    return pd.DataFrame(data)


def predict_ClinVar_order3(gpu_id=1, type='Insertion'):
    device = torch.device(f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu")
    transformer = train_model.load_model(device, model_dir='Model_Trained',
                                         model_name='pegRNA_Model_Merged_saved.order3_decoder_ori.pt')  # load model
    data_X = read_data_of_ClinVar(type)

    o = evaluate_model.transformer_predictor_order3(transformer, data_X, 1024, device)[0]

    pd.DataFrame({'Efficiency': o}).to_csv(f'Output/Efficiency.pegRNA.GRCh38.{type}.txt',
                                           sep='\t', index=False)


"""
@author hjs
"""


def read_data_of_Single(DF):
    """obtain the data in Article_293T for various types of regressions"""
    df = DF
    # df = pd.read_table(f'Data/pegRNA.GRCh38.{type}.txt', header=0)
    # df = df.loc[:, ['spacer sequence', 'RT', 'PBS', 'Original (Bold: sgRNA)', 'Efficiency(%)']]
    # df = df.iloc[:10, :]
    data = {'Target': [], 'Target_o2': [], 'Target_o3': [], 'RT': [], 'RT_o2': [], 'RT_o3': [],
            'PBS': [], 'PBS_o2': [], 'PBS_o3': []}

    len_RT = []
    len_PBS = []

    for i in range(len(df)):
        temp = df.iloc[i]
        print(temp)
        len_RT.append(len(temp['RT']))
        len_PBS.append(len(temp['PBS']))
    MAX_RT = max(len_RT)  # 31
    MAX_PBS = max(len_PBS)  # 17
    max_len_Target = 47
    print(f'Maximum length of (Target seq, PBS, RT): ({max_len_Target}, {MAX_PBS}, {MAX_RT})')

    id2char = list('ACGT')
    char2id = {char: i + 1 for i, char in enumerate(id2char)}
    char2id_o2 = {f'{char}{char_j}': i * len(id2char) + j + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char)}
    char2id_o3 = {f'{char}{char_j}{char_k}': i * len(id2char) * len(id2char) + j * len(id2char) + k + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char) for k, char_k in
                  enumerate(id2char)}
    for i, row in df.iterrows():
        seq = reverse_seq(complement_seq(row['PBS'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_PBS):
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_PBS - 1):
            temp.append(0)
        data['PBS_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_PBS - 2):
            temp.append(0)
        data['PBS_o3'].append(temp)

        seq = reverse_seq(complement_seq(row['RT'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_RT):
            temp.append(0)
        data['RT'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_RT - 1):
            temp.append(0)
        data['RT_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_RT - 2):
            temp.append(0)
        data['RT_o3'].append(temp)

        # seq = row['Original (Bold: sgRNA)'].upper()
        seq = row['Target(47bp)'].upper()
        # data['Target'].append([char2id[seq[j]] for j in range(0, len(seq))])
        # data['Target_o2'].append([char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)])
        # data['Target_o3'].append([char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)])
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_Target):
            temp.append(0)
        data['Target'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_Target - 1):
            temp.append(0)
        data['Target_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_Target - 2):
            temp.append(0)
        data['Target_o3'].append(temp)

        # data['Efficiency'].append(row['Efficiency(%)'] / 1)

    return pd.DataFrame(data)

def read_data_of_ClinVar_file(DF):
    """obtain the data in Article_293T for various types of regressions"""

    df = DF
    # df = df.loc[:, ['spacer sequence', 'RT', 'PBS', 'Original (Bold: sgRNA)', 'Efficiency(%)']]
    # df = df.iloc[:10, :]

    data = {'Target': [], 'Target_o2': [], 'Target_o3': [], 'RT': [], 'RT_o2': [], 'RT_o3': [],
            'PBS': [], 'PBS_o2': [], 'PBS_o3': []}

    len_RT = []
    len_PBS = []

    for i in range(len(df)):
        temp = df.iloc[i]
        print(temp['RT'],i)
        len_RT.append(len(temp['RT']))
        len_PBS.append(len(temp['PBS']))
    MAX_RT = max(len_RT)  # 31
    MAX_PBS = max(len_PBS)  # 17
    max_len_Target = 47
    print(f'Maximum length of (Target seq, PBS, RT): ({max_len_Target}, {MAX_PBS}, {MAX_RT})')

    id2char = list('ACGT')
    char2id = {char: i + 1 for i, char in enumerate(id2char)}
    char2id_o2 = {f'{char}{char_j}': i * len(id2char) + j + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char)}
    char2id_o3 = {f'{char}{char_j}{char_k}': i * len(id2char) * len(id2char) + j * len(id2char) + k + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char) for k, char_k in
                  enumerate(id2char)}
    for i, row in df.iterrows():
        seq = reverse_seq(complement_seq(row['PBS'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_PBS):
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_PBS - 1):
            temp.append(0)
        data['PBS_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_PBS - 2):
            temp.append(0)
        data['PBS_o3'].append(temp)

        seq = reverse_seq(complement_seq(row['RT'].upper()))
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), MAX_RT):
            temp.append(0)
        data['RT'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), MAX_RT - 1):
            temp.append(0)
        data['RT_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), MAX_RT - 2):
            temp.append(0)
        data['RT_o3'].append(temp)

        # seq = row['Original (Bold: sgRNA)'].upper()
        # seq = row['Target(47bp)'].upper()
        seq=row['Target'].upper()
        # data['Target'].append([char2id[seq[j]] for j in range(0, len(seq))])
        # data['Target_o2'].append([char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)])
        # data['Target_o3'].append([char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)])
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_Target):
            temp.append(0)
        data['Target'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_Target - 1):
            temp.append(0)
        data['Target_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_Target - 2):
            temp.append(0)
        data['Target_o3'].append(temp)

        # data['Efficiency'].append(row['Efficiency(%)'] / 1)

    return pd.DataFrame(data)


def optimal_combination():
    pass


def main():
    predict_ClinVar_order3(gpu_id=0, type='Insertion')


if __name__ == "__main__":
    main()
