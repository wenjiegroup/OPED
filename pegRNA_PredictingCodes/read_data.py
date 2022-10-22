#!/usr/bin/env python3

"""
author: Feng Liu
"""

import numpy as np
import pandas as pd
import re


one_hot_encoding = {'A': (1, 0, 0, 0),
                    'C': (0, 1, 0, 0),
                    'G': (0, 0, 1, 0),
                    'T': (0, 0, 0, 1)}


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


def read_data_for_sl():
    """obtain the data in NBT for various types of Shallow Learning methods"""

    data = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 1 (HT-training, test)',
                         header=1)
    # data = df.loc[:, :]

    max_len_PBS = max(data["PBS length"])  # 17
    max_len_RT = max(data["RT length"])  # 20
    max_len_PBS_RT = max(data["PBS-RT length"])  # 37
    max_len_Target = 47
    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    data_x = []
    data_y = []
    for i, row in data.iterrows():
        x_PBS = []
        temp = row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper()
        # temp = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][:row["PBS length"]])).upper()  # reversed and complementary sequence
        for j in range(len(temp)):
            x_PBS += one_hot_encoding[temp[j]]
        x_PBS += (0, 0, 0, 0) * (max_len_PBS - len(temp))
        x_RT = []
        temp = row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper()
        # temp = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][row["PBS length"]:])).upper()   # reversed and complementary sequence
        for j in range(len(temp)):
            x_RT += one_hot_encoding[temp[j]]
        x_RT += (0, 0, 0, 0) * (max_len_RT - len(temp))
        x_Target = []
        temp = row.iloc[2].upper()
        for j in range(max_len_Target):
            x_Target += one_hot_encoding[temp[j]]
        x_other = list(row[list(range(5, 8)) + list(range(9, 26))])
        x = x_RT + x_PBS + x_Target + x_other
        data_x.append(x)
        data_y.append(row['Measured PE efficiency'] / 1)
        # data_y.append(row['Measured PE efficiency'] if row['Measured PE efficiency'] >= 0 else 0)

    return np.array(data_x), np.array(data_y)


def read_data_for_sl_position_and_type(flag='Position'):   #flag= Position or Type
    """obtain the data in NBT for various types of Shallow Learning methods"""

    data = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 2 (Position, Type)',
                         header=1)
    data = data.replace('na', 0)
    # data = df.loc[:, :]

    max_len_PBS = max(data["PBS length"])  # 13
    max_len_RT = max(data["RT length"])  # 24
    max_len_PBS_RT = max(data["PBS-RT length"])  # 37
    max_len_Target = 47
    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    data_x = []
    data_y = []
    for i, row in data.iterrows():
        if not re.match(f'{flag}' + r'-\w+', row[0], re.I):
            continue
        x_PBS = []
        temp = row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper()
        # temp = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][:row["PBS length"]])).upper()  # reversed and complementary sequence
        for j in range(len(temp)):
            x_PBS += one_hot_encoding[temp[j]]
        x_PBS += (0, 0, 0, 0) * (max_len_PBS - len(temp))
        x_RT = []
        temp = row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper()
        # temp = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][row["PBS length"]:])).upper()   # reversed and complementary sequence
        for j in range(len(temp)):
            x_RT += one_hot_encoding[temp[j]]
        x_RT += (0, 0, 0, 0) * (max_len_RT - len(temp))
        x_Target = []
        temp = row.iloc[2].upper()
        for j in range(max_len_Target):
            x_Target += one_hot_encoding[temp[j]]
        x_other = list(row[list(range(5, 25))])
        x = x_RT + x_PBS + x_Target + x_other
        data_x.append(x)
        data_y.append(row['Measured PE efficiency'] / 1)
        # data_y.append(row['Measured PE efficiency'] if row['Measured PE efficiency'] >= 0 else 0)

    return np.array(data_x), np.array(data_y)


def read_data_for_rnn():
    """obtain the data in NBT for GRU"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 1 (HT-training, test)',
                       header=1)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(df["PBS length"])  # 17
    max_len_RT = max(df["RT length"])  # 20
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37
    max_len_Target = 47
    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    for i, row in df.iterrows():
        temp = [one_hot_encoding[s] for s in list(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())]
        # temp = [one_hot_encoding[s] for s in list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())))]  # reversed and complementary sequence of PBS
        for i in range(len(temp), max_len_PBS):
            temp.append((0, 0, 0, 0))
        data['PBS'].append(temp)
        temp = [one_hot_encoding[s] for s in list(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())]
        # temp = [one_hot_encoding[s] for s in list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())))]  # reversed and complementary sequence of RT
        for i in range(len(temp), max_len_RT):
            temp.append((0, 0, 0, 0))
        data['RT'].append(temp)
        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        data['Other'].append(tuple(row.iloc[list(range(5, 8)) + list(range(9, 26))]))
        data['Target'].append([one_hot_encoding[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)


def read_data_for_rnn_position_and_type(flag='Position'):
    """obtain the data in NBT for GRU"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 2 (Position, Type)',
                       header=1)
    df = df.replace('na', 0)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(df["PBS length"])  # 17
    max_len_RT = max(df["RT length"])  # 20
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37
    max_len_Target = 47
    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    for i, row in df.iterrows():
        if not re.match(f'{flag}' + r'-\w+', row[0], re.I):
            continue

        temp = [one_hot_encoding[s] for s in list(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())]
        # temp = [one_hot_encoding[s] for s in list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())))]  # reversed and complementary sequence of PBS
        for i in range(len(temp), max_len_PBS):
            temp.append((0, 0, 0, 0))
        data['PBS'].append(temp)
        temp = [one_hot_encoding[s] for s in list(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())]
        # temp = [one_hot_encoding[s] for s in list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())))]  # reversed and complementary sequence of RT
        for i in range(len(temp), max_len_RT):
            temp.append((0, 0, 0, 0))
        data['RT'].append(temp)
        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        data['Other'].append(tuple(row.iloc[list(range(5, 25))]))
        data['Target'].append([one_hot_encoding[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)


def read_data_of_for_transformer(max_len_Target=47, max_len_PBS=17, max_len_RT=20):
    """obtain the data in NBT for transformer"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 1 (HT-training, test)',
                       header=1)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'RT': [], 'PBS': [], 'Efficiency': []}
    # data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(max_len_PBS, max(df["PBS length"]))   # 17
    max_len_RT = max(max_len_RT, max(df["RT length"]))  # 20
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37

    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    id2char = list('ACGT')
    char2id = {char: i+1 for i, char in enumerate(id2char)}

    for i, row in df.iterrows():
        temp = [char2id[s] for s in list(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())]
        for j in range(len(temp), max_len_PBS):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id[s] for s in list(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())]
        for j in range(len(temp), max_len_RT):
            temp.append(0)
        data['RT'].append(temp)
        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        # data['Other'].append(list(row.iloc[list(range(5, 8)) + list(range(9, 26))]))
        data['Target'].append([char2id[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)


def read_data_for_transformer_position_and_type(flag='Position', max_len_Target=47, max_len_PBS=17, max_len_RT=20):   #flag= Position or Type
    """obtain the data in NBT for transformer"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 2 (Position, Type)',
                       header=1)
    df = df.replace('na', 0)
    # df.iloc[:, 10] = df.iloc[:, 10].replace('na', 48)
    # df.iloc[:, 12] = df.iloc[:, 12].replace('na', -7.2)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'RT': [], 'PBS': [], 'Efficiency': []}
    # data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(max_len_PBS, max(df["PBS length"]))  # max(df["PBS length"]) == 13, but 17 better than 13
    max_len_RT = max(max_len_RT, max(df["RT length"]))  # 24, max(df["RT length"]) == 24
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37

    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    id2char = list('ACGT')
    char2id = {char: i+1 for i, char in enumerate(id2char)}

    for i, row in df.iterrows():
        if not re.match(f'{flag}' + r'-\w+', row[0], re.I):
            continue
        # temp = [char2id[s] for s in list(row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper())]
        temp = [char2id[s] for s in
                list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'].upper()))[:row["PBS length"]])]
        for j in range(len(temp), max_len_PBS):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS'].append(temp)
        # temp = [char2id[s] for s in list(row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper())]
        temp = [char2id[s] for s in
                list(reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'].upper()))[row["PBS length"]:])]
        for j in range(len(temp), max_len_RT):
            temp.append(0)
        data['RT'].append(temp)
        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        # data['Other'].append(list(row.iloc[list(range(5, 25))]))
        data['Target'].append([char2id[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)


def read_data_of_for_transformer_order3(max_len_Target=47, max_len_PBS=17, max_len_RT=20):
    """obtain the data in NBT for transformer"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 1 (HT-training, test)',
                       header=1)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'Target_o2': [], 'Target_o3': [], 'RT': [], 'RT_o2': [], 'RT_o3': [],
            'PBS': [], 'PBS_o2': [], 'PBS_o3': [], 'Efficiency': []}
    # data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(max_len_PBS, max(df["PBS length"]))   # 17
    max_len_RT = max(max_len_RT, max(df["RT length"]))  # 20
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37

    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    id2char = list('ACGT')
    char2id = {char: i+1 for i, char in enumerate(id2char)}
    char2id_o2 = {f'{char}{char_j}': i * len(id2char) + j + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char)}
    char2id_o3 = {f'{char}{char_j}{char_k}': i * len(id2char) * len(id2char) + j * len(id2char) + k + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char) for k, char_k in enumerate(id2char)}

    for i, row in df.iterrows():
        seq = row['3\' extension sequence of pegRNA'][:row["PBS length"]].upper()
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_PBS):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id_o2[seq[j:j+2]] for j in range(0, len(seq)-1)]
        for j in range(len(temp), max_len_PBS-1):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_PBS - 2):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS_o3'].append(temp)

        seq = row['3\' extension sequence of pegRNA'][row["PBS length"]:].upper()
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_RT):
            temp.append(0)
        data['RT'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_RT - 1):
            # temp.insert(0, 0)
            temp.append(0)
        data['RT_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_RT - 2):
            # temp.insert(0, 0)
            temp.append(0)
        data['RT_o3'].append(temp)

        seq = row.iloc[2].upper()
        data['Target'].append([char2id[seq[j]] for j in range(0, len(seq))])
        data['Target_o2'].append([char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)])
        data['Target_o3'].append([char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)])

        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        # data['Other'].append(list(row.iloc[list(range(5, 8)) + list(range(9, 26))]))
        # data['Target'].append([char2id[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)


def read_data_for_transformer_position_and_type_order3(flag='Position', max_len_Target=47, max_len_PBS=17, max_len_RT=20):   #flag= Position or Type
    """obtain the data in NBT for transformer"""

    df = pd.read_excel('../Supplementary Table 4.xlsx', sheet_name='Library 2 (Position, Type)', header=1)
    df = df.replace('na', 0)
    # df.iloc[:, 10] = df.iloc[:, 10].replace('na', 48)
    # df.iloc[:, 12] = df.iloc[:, 12].replace('na', -7.2)
    # raw_data = df.iloc[:, [2, 4, 5, 26]]
    data = {'Target': [], 'Target_o2': [], 'Target_o3': [], 'RT': [], 'RT_o2': [], 'RT_o3': [],
            'PBS': [], 'PBS_o2': [], 'PBS_o3': [], 'Efficiency': []}
    # data = {'Target': [], 'RT': [], 'PBS': [], 'Efficiency': []}
    # data = {'Target': [], 'RT': [], 'PBS': [], 'Other': [], 'Efficiency': []}

    max_len_PBS = max(max_len_PBS, max(df["PBS length"]))  # max(df["PBS length"]) == 13, but 17 better than 13
    max_len_RT = max(max_len_RT, max(df["RT length"]))  # 24, max(df["RT length"]) == 24
    max_len_PBS_RT = max(df["PBS-RT length"])  # 37

    print(f'Maximum length of (Target, PBS, RT, PBS+RT): ({max_len_Target}, {max_len_PBS}, {max_len_RT}, {max_len_PBS_RT})')

    id2char = list('ACGT')
    char2id = {char: i + 1 for i, char in enumerate(id2char)}
    char2id_o2 = {f'{char}{char_j}': i * len(id2char) + j + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char)}
    char2id_o3 = {f'{char}{char_j}{char_k}': i * len(id2char) * len(id2char) + j * len(id2char) + k + 1
                  for i, char in enumerate(id2char) for j, char_j in enumerate(id2char) for k, char_k in
                  enumerate(id2char)}

    for i, row in df.iterrows():
        if not re.match(f'{flag}' + r'-\w+', row[0], re.I):
            continue

        seq = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'].upper()))[:row["PBS length"]]
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_PBS):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_PBS - 1):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_PBS - 2):
            # temp.insert(0, 0)
            temp.append(0)
        data['PBS_o3'].append(temp)

        seq = reverse_seq(complement_seq(row['3\' extension sequence of pegRNA'].upper()))[row["PBS length"]:]
        temp = [char2id[seq[j]] for j in range(0, len(seq))]
        for j in range(len(temp), max_len_RT):
            temp.append(0)
        data['RT'].append(temp)
        temp = [char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)]
        for j in range(len(temp), max_len_RT - 1):
            # temp.insert(0, 0)
            temp.append(0)
        data['RT_o2'].append(temp)
        temp = [char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)]
        for j in range(len(temp), max_len_RT - 2):
            # temp.insert(0, 0)
            temp.append(0)
        data['RT_o3'].append(temp)

        seq = row.iloc[2].upper()
        data['Target'].append([char2id[seq[j]] for j in range(0, len(seq))])
        data['Target_o2'].append([char2id_o2[seq[j:j + 2]] for j in range(0, len(seq) - 1)])
        data['Target_o3'].append([char2id_o3[seq[j:j + 3]] for j in range(0, len(seq) - 2)])

        data['Efficiency'].append(row['Measured PE efficiency'] / 1)
        # data['Other'].append(list(row.iloc[list(range(5, 25))]))
        # data['Target'].append([char2id[s] for s in list(row.iloc[2].upper())])

    return pd.DataFrame(data)

