#!/usr/bin/env python3

"""
author: Feng Liu
"""

import numpy as np
import pandas as pd
import time
import math
import random
import sys
import os
import copy

import torch
import torch.nn as nn
import torch.nn.functional as F

from pegRNA_PredictingCodes import read_data,evaluate_model
# from pegsub import *
#wenjianmulujiegou hunluan
sys.path.insert(0,'./pegRNA_PredictingCodes')  #huifujiegou

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)  # for CPU
    torch.cuda.manual_seed(seed)  # for this GPU
    torch.cuda.manual_seed_all(seed)  # for all GPU


def resample_with_replacement(X_train, y_train, n):
    # resample with replace to the maximum number

    print(f'resample with replace to #{n}')
    # ind = np.array(range(X_train.shape[0]))
    # ind = np.concatenate((np.random.choice(ind, n), ind), 0)
    ind = np.random.choice(X_train.shape[0], n)
    X_train = X_train.iloc[ind, :]
    y_train = y_train.iloc[ind]

    return X_train, y_train


def train_and_test_sl(X_train, X_test, y_train, y_test):
    """train and test various types of Shallow Learning methods"""

    from sklearn.svm import NuSVR
    from sklearn.neural_network import MLPRegressor
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.ensemble import AdaBoostRegressor
    from sklearn.ensemble import GradientBoostingRegressor
    from xgboost import XGBRegressor

    r_of_model = {'pearson': {}, 'spearman': {}}
    start = time.time()

    # model 1: Support Vector Regression
    m1 = NuSVR()
    m1.fit(X_train, y_train)
    m_name = "SVM"
    temp = evaluate_model.evaluate_sl(m1, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t1 = time.time()
    print(f'{m_name} time costed: {t1 - start}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 2: Multi-layer Perceptron Regression
    m2 = MLPRegressor(hidden_layer_sizes=(500,), activation='logistic', solver='adam',
                      learning_rate='adaptive', learning_rate_init=0.001, max_iter=200)
    m2.fit(X_train, y_train)
    m_name = "MLP"
    temp = evaluate_model.evaluate_sl(m2, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t2 = time.time()
    print(f'{m_name} time costed: {t2 - t1}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 3: K-Nearest Neighbors Regression
    m3 = KNeighborsRegressor(n_neighbors=5)
    m3.fit(X_train, y_train)
    m_name = "KNN"
    temp = evaluate_model.evaluate_sl(m3, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t3 = time.time()
    print(f'{m_name} time costed: {t3 - t2}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 4: Decision Tree Regression
    m4 = DecisionTreeRegressor()
    m4.fit(X_train, y_train)
    m_name = "DT"
    temp = evaluate_model.evaluate_sl(m4, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t4 = time.time()
    print(f'{m_name} time costed: {t4 - t3}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 5: Random Forest Regression
    m5 = RandomForestRegressor()
    m5.fit(X_train, y_train)
    m_name = "RF"
    temp = evaluate_model.evaluate_sl(m5, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t5 = time.time()
    print(f'{m_name} time costed: {t5 - t4}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 6: AdaBoost Regression
    m6 = AdaBoostRegressor(DecisionTreeRegressor(), n_estimators=50)
    m6.fit(X_train, y_train)
    m_name = "AdaBoost"
    temp = evaluate_model.evaluate_sl(m6, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t6 = time.time()
    print(f'{m_name} time costed: {t6 - t5}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 7: Gradient Boosted Decision Trees Regression
    m7 = GradientBoostingRegressor()
    m7.fit(X_train, y_train)
    m_name = "GBDT"
    temp = evaluate_model.evaluate_sl(m7, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t7 = time.time()
    print(f'{m_name} time costed: {t7 - t6}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    # model 8: XGBoost Regression
    m8 = XGBRegressor(objective='reg:squarederror')
    m8.fit(X_train, y_train)
    m_name = "XGBoost"
    temp = evaluate_model.evaluate_sl(m8, X_train, X_test, y_train, y_test)
    r_of_model['pearson'][m_name] = temp['pearson']
    r_of_model['spearman'][m_name] = temp['spearman']
    t8 = time.time()
    print(f'{m_name} time costed: {t8 - t7}s')
    # print(m_name, r_of_model['pearson'][m_name], r_of_model['spearman'][m_name])

    end = time.time()
    print(f'total time costed: {end - start}s')
    return {'pearson': pd.DataFrame(r_of_model['pearson'], index=('train', 'test')),
            'spearman': pd.DataFrame(r_of_model['spearman'], index=('train', 'test'))}


class EncoderRNN(nn.Module):
    def __init__(self, input_size, hidden_size, hidden_size_fully, output_size, drop_out_p, bidirectional=True,
                 num_layers=1):
        """
        :param input_size: a tuple of (input_size_original, input_size_pbs, input_size_rt)
        :param hidden_size: a tuple of (hidden_size_original, hidden_size_pbs, hidden_size_rt)
        """
        super(EncoderRNN, self).__init__()
        self.hidden_size = hidden_size
        self.bidirectional = bidirectional
        self.num_layers = num_layers
        self.gru_original = nn.GRU(input_size[0], hidden_size[0], num_layers=num_layers,
                                   batch_first=True, bidirectional=bidirectional)
        self.gru_pbs = nn.GRU(input_size[1], hidden_size[1], num_layers=num_layers,
                              batch_first=True, bidirectional=bidirectional)
        self.gru_rt = nn.GRU(input_size[2], hidden_size[2], num_layers=num_layers,
                             batch_first=True, bidirectional=bidirectional)
        self.linear_att_original = nn.Linear(hidden_size[0] * (2 if bidirectional else 1), 1)
        self.linear_att_pbs = nn.Linear(hidden_size[1] * (2 if bidirectional else 1), 1)
        self.linear_att_rt = nn.Linear(hidden_size[2] * (2 if bidirectional else 1), 1)
        self.fully_connected_layers = nn.Sequential(
            nn.Linear(sum(hidden_size) * (2 if bidirectional else 1) + input_size[3], hidden_size_fully[0]),
            nn.ReLU(),
            nn.Dropout(drop_out_p),
            nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
            nn.ReLU(),
            nn.Dropout(drop_out_p),
            nn.Linear(hidden_size_fully[1], output_size),
        )
        # self.fully_connected_layers = nn.Sequential(
        # nn.Linear(sum(hidden_size) * (2 if bidirectional else 1), hidden_size_fully[0]),
        # nn.ReLU(),
        # nn.Dropout(drop_out_p),
        # nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
        # nn.ReLU(),
        # nn.Dropout(drop_out_p),
        # nn.Linear(hidden_size_fully[1], output_size),
        # nn.Dropout(drop_out_p)
        # )
        # self.linear_combined = nn.Linear(sum(hidden_size) * (2 if bidirectional else 1), output_size)
        # self.softmax = nn.Softmax(dim=1)

    def forward(self, input, hidden):
        """
        :param input: a tuple of (input_original, input_pbs, input_rt)
        :param hidden:
        :return:
        """

        # original sequence
        output_original = input[0]
        hidden_original = hidden[0]
        output_original, hidden_original = self.gru_original(output_original, hidden_original)
        scores_original = self.linear_att_original(output_original).squeeze(-1)
        scores_original = F.softmax(scores_original, dim=1)
        output_original = scores_original.unsqueeze(-2).bmm(output_original).squeeze(-2)

        # pbs
        output_pbs = input[1]
        hidden_pbs = hidden[1]
        output_pbs, hidden_pbs = self.gru_pbs(output_pbs, hidden_pbs)
        scores_pbs = self.linear_att_pbs(output_pbs).squeeze(-1)
        scores_pbs = F.softmax(scores_pbs, dim=1)
        output_pbs = scores_pbs.unsqueeze(-2).bmm(output_pbs).squeeze(-2)

        # rt
        output_rt = input[2]
        hidden_rt = hidden[2]
        output_rt, hidden_rt = self.gru_rt(output_rt, hidden_rt)
        scores_rt = self.linear_att_rt(output_rt).squeeze(-1)
        scores_rt = F.softmax(scores_rt, dim=1)
        output_rt = scores_rt.unsqueeze(-2).bmm(output_rt).squeeze(-2)

        # other 20 features
        other = input[3]

        # combined = torch.cat((output_original[:, -1, :], output_pbs[:, -1, :], output_rt[:, -1, :], other), 1)
        combined = torch.cat((output_original, output_pbs, output_rt, other), 1)
        output = (self.fully_connected_layers(combined))
        # output = (self.linear_combined(combined))

        return output

    def initHidden(self, bacth_size, device):
        return torch.zeros(self.num_layers * (2 if self.bidirectional else 1),
                           bacth_size, self.hidden_size[0], device=device)


class PositionalEncoding(nn.Module):
    r"""PositionEncoding without parameters is faster, and smaller.
        It's said that similar to it with parameters in some tasks.
        Inject some information about the relative or absolute position of the tokens
        in the sequence. The positional encodings have the same dimension as
        the embeddings, so that the two can be summed. Here, we use sine and cosine
        functions of different frequencies.
    .. math::
        \text{PosEncoder}(pos, 2i) = sin(pos/10000^(2i/d_model))
        \text{PosEncoder}(pos, 2i+1) = cos(pos/10000^(2i/d_model))
        \text{where pos is the word position and i is the embed idx)
    Args:
        d_model: the embed dim (required).
        dropout: the dropout value (default=0.1).
        max_len: the max. length of the incoming sequence (default=5000).
    Examples:
        >>> pos_encoder = PositionalEncoding(d_model)
    """

    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)

        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)

    def forward(self, x):
        r"""Inputs of forward function
        Args:
            x: the sequence fed to the positional encoder model (required).
        Shape:
            x: [sequence length, batch size, embed dim]
            output: [sequence length, batch size, embed dim]
        Examples:
            >>> output = pos_encoder(x)
        """

        x = x + self.pe[:x.size(0), :]
        return self.dropout(x)


class LearnedPositionEncoding(nn.Embedding):
    """for some tasks, PositionEncoding with parameters are clearly better that one without.
    PositionEncoding layer with parameters is easy to definer, inherit from a nn.Embedding, add a dropout
    because nn.Embedding contain A weight matrix that can be oriented by index
    """
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super().__init__(max_len, d_model)
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, x):
        weight = self.weight.data.unsqueeze(1)
        x = x + weight[:x.size(0), :]
        return self.dropout(x)


class TransformerEncoderModel(nn.Module):
    """Container module with an encoder, a recurrent or transformer module, and a decoder."""

    def __init__(self, ntoken=4, embedding_size=512, hidden_size=[2048, 2048, 2048], hidden_size_fully=None,
                 output_size=1, nhead=8, num_encoder_layers=[6, 6, 6], dropout=0.1, softmax_bool=False,
                 other_size=0):
        super(TransformerEncoderModel, self).__init__()
        self.model_type = 'Transformer'
        self.embedding_size = embedding_size
        self.nhead = nhead
        self.n = 3
        self.other_size = other_size
        self.dropout = dropout

        self.embedding = nn.Embedding(ntoken+1, embedding_size, padding_idx=0)
        self.pos_encoder = PositionalEncoding(embedding_size, dropout)

        self.softmax_bool = softmax_bool
        self.encoder = nn.ModuleList()
        self.linear_att = nn.ModuleList()

        for i in range(0, self.n):
            encoder_i, linear_att_i = self.init_transformer_encoder(d_model=embedding_size, nhead=nhead,
                                                                    dim_feedforward=hidden_size[i],
                                                                    num_encoder_layer=num_encoder_layers[i],
                                                                    dropout=dropout)
            self.encoder.append(encoder_i)
            self.linear_att.append(linear_att_i)
        if hidden_size_fully is None:
            self.fully_connected_layers = nn.Linear(self.n * embedding_size + self.other_size, output_size)
        elif len(hidden_size_fully) == 1:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n * embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], output_size), )
        elif len(hidden_size_fully) == 2:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n * embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[0]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[1]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], output_size), )
        else:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n * embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], hidden_size_fully[2]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[2], output_size), )
        # self.init_weights()
        self._reset_parameters()

    def init_transformer_encoder(self, d_model, nhead, dim_feedforward, num_encoder_layer, dropout):
        try:
            from torch.nn import Transformer, TransformerEncoder, TransformerEncoderLayer, LayerNorm
        except Exception:
            raise ImportError('TransformerEncoder module does not exist in PyTorch 1.1 or lower.')
        encoder_layer = TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward, dropout=dropout)
        encoder_norm = LayerNorm(d_model)
        encoder = TransformerEncoder(encoder_layer, num_encoder_layer, encoder_norm)
        # encoder = TransformerEncoder(encoder_layer, num_encoder_layer, norm=None)
        linear_att = nn.Linear(self.embedding_size, 1)
        # linear_att = nn.Linear(self.embedding_size, self.embedding_size)
        return encoder, linear_att

    def _reset_parameters(self):
        r"""Initiate parameters in the transformer model."""

        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    # def init_weights(self):
    #     initrange = 0.1
    #     nn.init.uniform_(self.embedding.weight, -initrange, initrange)
    #     nn.init.uniform_(self.fully_connected_layers.weight, -initrange, initrange)

    def padding_mask(self, seq, pad_idx):
        return seq == pad_idx   #.unsqueeze(-2)

    def forward(self, input):

        memories = []
        att_weight = []
        for i in range(0, len(self.encoder)):
            seq = self.embedding(input[i].T) * math.sqrt(self.embedding_size)
            seq = self.pos_encoder(seq)
            mask = self.padding_mask(input[i], pad_idx=0)
            # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask))
            # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask).masked_fill(mask.T.unsqueeze(-1), 0))
            values = self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask) #.masked_fill(mask.T.unsqueeze(-1), 0)
            # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) / math.sqrt(self.embedding_size)  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2) / math.sqrt(self.embedding_size)  # values * mask.T
            # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) * 1e6  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2)  # values * mask.T
            # scores = torch.mean(self.linear_att[i](values) * values, 2)
            scores = self.linear_att[i](values).squeeze(-1)
            scores.masked_fill_(mask.T, float('-inf'))
            scores = F.softmax(scores, dim=0)
            att_weight.append(scores.T)
            memories.append(scores.unsqueeze(-1).permute(1, 2, 0).bmm(values.permute(1, 0, 2)).squeeze(-2)) # N*1*S @ N*S*E = N*1*E

        # combined = torch.cat((memories[0][0], memories[1][0], memories[2][0]), 1)
        # combined = torch.cat((torch.mean(memories[0], 0), torch.mean(memories[1], 0), torch.mean(memories[2], 0)), 1)
        combined = torch.cat((memories[0], memories[1], memories[2]), 1)
        # combined = torch.cat((memories[0], memories[1], memories[2], input[3]), 1)
        output = (self.fully_connected_layers(combined))
        if self.softmax_bool:
            output = F.softmax(output, dim=-1)

        return output, att_weight


class TransformerEncoderModelOrder3(nn.Module):
    """Container module with an encoder, a recurrent or transformer module, and a decoder."""

    def __init__(self, ntokens=[4, 16, 64], embedding_size=512, hidden_size=[2048, 2048, 2048], hidden_size_fully=None,
                 output_size=1, nhead=8, num_encoder_layers=[6, 6, 6], dropout=0.1, softmax_bool=False,
                 other_size=0):
        super(TransformerEncoderModelOrder3, self).__init__()
        self.model_type = 'Transformer'
        self.embedding_size = embedding_size
        self.nhead = nhead
        self.n = 3  # target/PBS/RT
        self.ntokens = ntokens  # order 1/2/3
        self.other_size = other_size
        self.dropout = dropout

        self.embedding = nn.ModuleList()
        self.pos_encoder = PositionalEncoding(embedding_size, dropout)

        self.softmax_bool = softmax_bool
        self.encoder = nn.ModuleList()
        self.linear_att = nn.ModuleList()

        for j in range(0, len(self.ntokens)):
            self.embedding.append(nn.Embedding(ntokens[j] + 1, embedding_size, padding_idx=0))
            for i in range(0, self.n):
                encoder_i, linear_att_i = self.init_transformer_encoder(d_model=embedding_size, nhead=nhead,
                                                                        dim_feedforward=hidden_size[i],
                                                                        num_encoder_layer=num_encoder_layers[i],
                                                                        dropout=dropout)
                self.encoder.append(encoder_i)
                self.linear_att.append(linear_att_i)
        if hidden_size_fully is None:
            self.fully_connected_layers = nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, output_size)
        elif len(hidden_size_fully) == 1:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], output_size), )
        elif len(hidden_size_fully) == 2:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[0]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[1]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], output_size), )
        else:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], hidden_size_fully[2]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[2], output_size), )
        # self.init_weights()
        self._reset_parameters()

    def init_transformer_encoder(self, d_model, nhead, dim_feedforward, num_encoder_layer, dropout):
        try:
            from torch.nn import Transformer, TransformerEncoder, TransformerEncoderLayer, LayerNorm
        except Exception:
            raise ImportError('TransformerEncoder module does not exist in PyTorch 1.1 or lower.')
        encoder_layer = TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward, dropout=dropout)
        encoder_norm = LayerNorm(d_model)
        encoder = TransformerEncoder(encoder_layer, num_encoder_layer, encoder_norm)
        # encoder = TransformerEncoder(encoder_layer, num_encoder_layer, norm=None)
        linear_att = nn.Linear(self.embedding_size, 1)
        # linear_att = nn.Linear(self.embedding_size, self.embedding_size)
        return encoder, linear_att

    def _reset_parameters(self):
        r"""Initiate parameters in the transformer model."""

        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    # def init_weights(self):
    #     initrange = 0.1
    #     nn.init.uniform_(self.embedding.weight, -initrange, initrange)
    #     nn.init.uniform_(self.fully_connected_layers.weight, -initrange, initrange)

    def padding_mask(self, seq, pad_idx):
        return seq == pad_idx   #.unsqueeze(-2)

    def forward(self, input):

        memories = []
        att_weight = []
        for j in range(0, len(self.ntokens)):
            for i in range(0, self.n):
                seq = self.embedding[j](input[j*len(self.ntokens)+i].T) * math.sqrt(self.embedding_size)
                seq = self.pos_encoder(seq)
                mask = self.padding_mask(input[j*len(self.ntokens)+i], pad_idx=0)
                # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask))
                # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask).masked_fill(mask.T.unsqueeze(-1), 0))
                values = self.encoder[j*len(self.ntokens)+i](src=seq, mask=None, src_key_padding_mask=mask) #.masked_fill(mask.T.unsqueeze(-1), 0)
                # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) / math.sqrt(self.embedding_size)  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2) / math.sqrt(self.embedding_size)  # values * mask.T
                # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) * 1e6  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2)  # values * mask.T
                # scores = torch.mean(self.linear_att[i](values) * values, 2)
                scores = self.linear_att[j*len(self.ntokens)+i](values).squeeze(-1)
                scores.masked_fill_(mask.T, float('-inf'))
                scores = F.softmax(scores, dim=0)
                att_weight.append(scores.T)
                memories.append(scores.unsqueeze(-1).permute(1, 2, 0).bmm(values.permute(1, 0, 2)).squeeze(-2)) # N*1*S @ N*S*E = N*1*E

        # combined = torch.cat((memories[0][0], memories[1][0], memories[2][0]), 1)
        # combined = torch.cat((torch.mean(memories[0], 0), torch.mean(memories[1], 0), torch.mean(memories[2], 0)), 1)
        combined = torch.cat([memories[i] for i in range(0, len(memories))], 1)
        # combined = torch.cat((memories[0], memories[1], memories[2], input[3]), 1)
        output = (self.fully_connected_layers(combined))
        if self.softmax_bool:
            output = F.softmax(output, dim=-1)

        return output, att_weight


class TransformerEncoderDecoderModelOrder3(nn.Module):
    """Container module with an encoder, a recurrent or transformer module, and a decoder."""

    def __init__(self, ntokens=[4, 16, 64], embedding_size=512, hidden_size=[2048, 2048, 2048], hidden_size_fully=None,
                 output_size=1, nhead=8, num_encoder_layers=[6, 6, 6], dropout=0.1, softmax_bool=False,
                 other_size=0):
        super(TransformerEncoderDecoderModelOrder3, self).__init__()
        self.model_type = 'TransformerEncoderDecoder'
        self.embedding_size = embedding_size
        self.nhead = nhead
        self.n = 3  # target/PBS/RT
        self.ntokens = ntokens  # order 1/2/3
        self.other_size = other_size
        self.dropout = dropout

        self.embedding = nn.ModuleList()
        self.pos_encoder = PositionalEncoding(embedding_size, dropout)

        self.softmax_bool = softmax_bool
        self.encoder_decoder = nn.ModuleList()
        self.linear_att = nn.ModuleList()
        print(self.model_type)

        for j in range(0, len(self.ntokens)):
            self.embedding.append(nn.Embedding(ntokens[j] + 1, embedding_size, padding_idx=0))
            for i in range(0, self.n):
                if i == 0:  # target
                    encoder_i, linear_att_i = self.init_transformer_encoder(d_model=embedding_size, nhead=nhead,
                                                                            dim_feedforward=hidden_size[i],
                                                                            num_encoder_layer=num_encoder_layers[i],
                                                                            dropout=dropout)
                else:   # PBS/RT
                    encoder_i, linear_att_i = self.init_transformer_decoder(d_model=embedding_size, nhead=nhead,
                                                                            dim_feedforward=hidden_size[i],
                                                                            num_decoder_layer=num_encoder_layers[i],
                                                                            dropout=dropout)
                self.encoder_decoder.append(encoder_i)
                self.linear_att.append(linear_att_i)
        if hidden_size_fully is None:
            self.fully_connected_layers = nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, output_size)
        elif len(hidden_size_fully) == 1:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], output_size), )
        elif len(hidden_size_fully) == 2:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[0]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[1]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], output_size), )
        elif len(hidden_size_fully) == 3:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n*len(self.ntokens)*embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], hidden_size_fully[2]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[2], output_size), )
        elif len(hidden_size_fully) == 4:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n * len(self.ntokens) * embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], hidden_size_fully[2]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[2], hidden_size_fully[3]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[3], output_size), )
        else:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(self.n * len(self.ntokens) * embedding_size + self.other_size, hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], hidden_size_fully[2]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[2], hidden_size_fully[3]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[3], hidden_size_fully[4]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[4], output_size), )
        # self.init_weights()
        self._reset_parameters()

    def init_transformer_encoder(self, d_model, nhead, dim_feedforward, num_encoder_layer, dropout):
        try:
            from torch.nn import Transformer, TransformerEncoder, TransformerEncoderLayer, LayerNorm
        except Exception:
            raise ImportError('TransformerEncoder module does not exist in PyTorch 1.1 or lower.')
        encoder_layer = TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward, dropout=dropout)
        encoder_norm = LayerNorm(d_model)
        encoder = TransformerEncoder(encoder_layer, num_encoder_layer, encoder_norm)
        # encoder = TransformerEncoder(encoder_layer, num_encoder_layer, norm=None)
        linear_att = nn.Linear(self.embedding_size, 1)
        # linear_att = nn.Linear(self.embedding_size, self.embedding_size)
        return encoder, linear_att

    def init_transformer_decoder(self, d_model, nhead, dim_feedforward, num_decoder_layer, dropout):
        try:
            from torch.nn import Transformer, TransformerDecoder, TransformerDecoderLayer, LayerNorm
        except Exception:
            raise ImportError('TransformerDecoder module does not exist in PyTorch 1.1 or lower.')
        decoder_layer = TransformerDecoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward, dropout=dropout)
        decoder_norm = LayerNorm(d_model)
        decoder = TransformerDecoder(decoder_layer, num_decoder_layer, decoder_norm)
        # decoder = TransformerDecoder(decoder_layer, num_decoder_layer, norm=None)
        linear_att = nn.Linear(self.embedding_size, 1)
        # linear_att = nn.Linear(self.embedding_size, self.embedding_size)
        return decoder, linear_att

    def _reset_parameters(self):
        r"""Initiate parameters in the transformer model."""

        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    # def init_weights(self):
    #     initrange = 0.1
    #     nn.init.uniform_(self.embedding.weight, -initrange, initrange)
    #     nn.init.uniform_(self.fully_connected_layers.weight, -initrange, initrange)

    def padding_mask(self, seq, pad_idx):
        return seq == pad_idx   #.unsqueeze(-2)

    def forward(self, input):

        memories = []
        att_weight = []
        encoder_memories = []
        encoder_masks = []
        for j in range(0, len(self.ntokens)):
            for i in range(0, self.n):
                seq = self.embedding[j](input[j*len(self.ntokens)+i].T) * math.sqrt(self.embedding_size)
                seq = self.pos_encoder(seq)
                mask = self.padding_mask(input[j*len(self.ntokens)+i], pad_idx=0)
                # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask))
                # memories.append(self.encoder[i](src=seq, mask=None, src_key_padding_mask=mask).masked_fill(mask.T.unsqueeze(-1), 0))
                if i == 0:  # target
                    values = self.encoder_decoder[j*len(self.ntokens)+i](src=seq, mask=None, src_key_padding_mask=mask) #.masked_fill(mask.T.unsqueeze(-1), 0)
                    encoder_memories.append(values)
                    encoder_masks.append(mask)
                else:   # PBS/RT
                    values = self.encoder_decoder[j * len(self.ntokens) + i](tgt=seq, memory=encoder_memories[-1],
                                                                             tgt_key_padding_mask=mask,
                                                                             memory_key_padding_mask=encoder_masks[-1])
                # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) / math.sqrt(self.embedding_size)  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2) / math.sqrt(self.embedding_size)  # values * mask.T
                # scores = torch.sum(values.masked_fill(mask.T.unsqueeze(-1), 0), 2) * 1e6  #same with torch.sum(values * ~mask.T.unsqueeze(-1), 2)  # values * mask.T
                # scores = torch.mean(self.linear_att[i](values) * values, 2)
                scores = self.linear_att[j*len(self.ntokens)+i](values).squeeze(-1)
                scores.masked_fill_(mask.T, float('-inf'))
                scores = F.softmax(scores, dim=0)
                att_weight.append(scores.T)
                memories.append(scores.unsqueeze(-1).permute(1, 2, 0).bmm(values.permute(1, 0, 2)).squeeze(-2)) # N*1*S @ N*S*E = N*1*E

        encoder_memories = []
        encoder_masks = []

        # combined = torch.cat((memories[0][0], memories[1][0], memories[2][0]), 1)
        # combined = torch.cat((torch.mean(memories[0], 0), torch.mean(memories[1], 0), torch.mean(memories[2], 0)), 1)
        combined = torch.cat([memories[i] for i in range(0, len(memories))], 1)
        # combined = torch.cat((memories[0], memories[1], memories[2], input[3]), 1)
        output = (self.fully_connected_layers(combined))
        if self.softmax_bool:
            output = F.softmax(output, dim=-1)

        return output, att_weight


def train_rnn(rnn, criterion, optimizer, X_train, y_train, epoch_num, batch_size, input_size, hidden_size,
              hidden_size_fully, num_layers, output_size, drop_out_p, weight_decay, bidirectional, device):
    # train
    rnn.train()
    n = len(X_train)
    print(f'epoch_num = {epoch_num}\nbatch_size = {batch_size}\ninput_size = {input_size}\nhidden_size = {hidden_size}'
          f'\nhidden_size_fully = {hidden_size_fully}\nnum_layers = {num_layers}\noutput_size = {output_size}'
          f'\ndrop_out_p = {drop_out_p}\nweight_decay = {weight_decay}\nbidirectional = {bidirectional}')
    print(f'Training data: {n}')
    batch_num = n // batch_size
    index_shuffled = np.random.permutation(n)
    X_train = X_train.iloc[index_shuffled, :]
    y_train = y_train.iloc[index_shuffled]
    # train_ds = TensorDataset(torch.tensor(np.array(X_train)), torch.tensor(np.array(y_train)))
    # train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    print('Start training...')
    start = time.time()
    for epoch in range(epoch_num):
        if epoch % 10 == 0:
            print(f'    Epoch {epoch}:')
        for i in range(batch_num):
            start_i = i * batch_size
            end_i = start_i + batch_size
            xb = X_train.iloc[start_i:end_i, :]
            yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.float32),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.float32),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.float32),
                     torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            hidden0 = (rnn.initHidden(len(xb['PBS']), device),) * 3

            optimizer.zero_grad()
            outputs = rnn(input, hidden0)
            loss = criterion(outputs.squeeze(), torch.tensor(list(yb), device=device, dtype=torch.float32))
            loss.backward()
            # torch.nn.utils.clip_grad_norm_(rnn.parameters(), 0.25)
            optimizer.step()
    print(f'Training time:{time.time() - start}')


def train_transformer(transformer, criterion, optimizer, scheduler, X_train, X_test, y_train, y_test,
                      epoch_num, batch_size, device, best_epoch=True):
    best_model_dict = copy.deepcopy(transformer.state_dict())
    best_r = 0.0
    best_e = 0

    # train
    n = len(X_train)
    print(f'Training data: {n}')
    batch_num = n // batch_size + 1
    index_shuffled = np.random.permutation(n)
    X_train = X_train.iloc[index_shuffled, :]
    y_train = y_train.iloc[index_shuffled]
    # train_ds = TensorDataset(torch.tensor(np.array(X_train)), torch.tensor(np.array(y_train)))
    # train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    print('Epoch:', end="")
    start = time.time()
    for epoch in range(epoch_num):
        transformer.train()
        if epoch % 10 == 0:
            # print(f'\t{epoch}', end="")
            sys.stdout.write(f'  {epoch}')
            sys.stdout.flush()
        for i in range(batch_num):
            start_i = i * batch_size
            end_i = start_i + batch_size
            xb = X_train.iloc[start_i:end_i, :]
            yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long))
                     # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))

            optimizer.zero_grad()
            outputs, _ = transformer(input)
            loss = criterion(outputs.squeeze(), torch.tensor(list(yb), device=device, dtype=torch.float32))
            loss.backward()
            # torch.nn.utils.clip_grad_norm_(rnn.parameters(), 0.25)
            optimizer.step()

        scheduler.step()
        epoch_r = evaluate_model.evaluate_transformer(transformer, X_test, y_test, 1024, device, False)[0]['pearson'][0]
        if epoch_r > best_r:
            best_model_dict = copy.deepcopy(transformer.state_dict())
            best_r = epoch_r
            best_e = epoch
    if best_epoch:
        transformer.load_state_dict(best_model_dict)
        print(f'\nBest epoch: {best_e}')
        print(f'Training time:{time.time() - start:.2f}s')
    else:
        print(f'\nTraining time:{time.time() - start:.2f}s')

    return transformer


def train_transformer_order3(transformer, criterion, optimizer, scheduler, X_train, X_test, y_train, y_test,
                             epoch_num, batch_size, device, best_epoch=True):
    best_model_dict = copy.deepcopy(transformer.state_dict())
    best_r = 0.0
    best_e = 0

    # train
    n = len(X_train)
    print(f'Training data: {n}')
    batch_num = n // batch_size + 1
    index_shuffled = np.random.permutation(n)
    X_train = X_train.iloc[index_shuffled, :]
    y_train = y_train.iloc[index_shuffled]
    # train_ds = TensorDataset(torch.tensor(np.array(X_train)), torch.tensor(np.array(y_train)))
    # train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    print('Epoch:', end="")
    start = time.time()
    for epoch in range(epoch_num):
        transformer.train()
        if epoch % 10 == 0:
            # print(f'\t{epoch}', end="")
            sys.stdout.write(f'  {epoch}')
            sys.stdout.flush()
        for i in range(batch_num):
            start_i = i * batch_size
            end_i = start_i + batch_size
            xb = X_train.iloc[start_i:end_i, :]
            yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long),)
                     # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))

            optimizer.zero_grad()
            outputs, _ = transformer(input)
            loss = criterion(outputs.squeeze(), torch.tensor(list(yb), device=device, dtype=torch.float32))
            loss.backward()
            # torch.nn.utils.clip_grad_norm_(rnn.parameters(), 0.25)
            optimizer.step()

        scheduler.step()
        epoch_r = evaluate_model.evaluate_transformer_order3(transformer, X_test, y_test, 1024, device, False)[0]['pearson'][0]
        if epoch_r > best_r:
            best_model_dict = copy.deepcopy(transformer.state_dict())
            best_r = epoch_r
            best_e = epoch
    if best_epoch:
        transformer.load_state_dict(best_model_dict)
        print(f'\nBest epoch: {best_e}')
        print(f'Training time:{time.time() - start:.2f}s')
    else:
        print(f'\nTraining time:{time.time() - start:.2f}s')

    return transformer


def train_and_test_rnn(X_train, X_test, y_train, y_test, device):
    epoch_num = 100 # 10 for lib1 and 100 for lib2
    batch_size = 128
    input_size = (4, 4, 4, 20)
    hidden_size = (50, 50, 50)
    hidden_size_fully = (100, 10)
    drop_out_p = 0
    num_layers = 2
    output_size = 1
    weight_decay = 0
    bidirectional = True
    rnn = EncoderRNN(input_size, hidden_size, hidden_size_fully, output_size, drop_out_p, bidirectional = bidirectional,
                     num_layers=num_layers).to(device)
    criterion = torch.nn.MSELoss()
    # optimizer = optim.SGD(rnn.parameters(), lr=0.5, momentum=0.9)
    # optimizer = torch.optim.RMSprop(rnn.parameters(), lr=0.01)
    optimizer = torch.optim.Adam(rnn.parameters(), lr=0.001, weight_decay=weight_decay)

    train_rnn(rnn, criterion, optimizer, X_train, y_train, epoch_num, batch_size, input_size, hidden_size,
              hidden_size_fully, num_layers, output_size, drop_out_p, weight_decay, bidirectional, device)
    coefficient = evaluate_model.evaluate_rnn(rnn, X_train, y_train, X_test, y_test)

    return coefficient, rnn


class CustomSchedule(torch.optim.lr_scheduler._LRScheduler):
    def __init__(self, optimizer, d_model, warm_steps=4):
        self.optimizer = optimizer
        self.d_model = d_model
        self.warmup_steps = warm_steps

        super(CustomSchedule, self).__init__(optimizer)

    def get_lr(self):
        """
        # rsqrt: 1 / sqrt{x}
        arg1 = torch.rsqrt(torch.tensor(self._step_count, dtype=torch.float32))
        arg2 = torch.tensor(self._step_count * (self.warmup_steps ** -1.5), dtype=torch.float32)
        dynamic_lr = torch.rsqrt(self.d_model) * torch.minimum(arg1, arg2)
        """
        # print('*'*27, self._step_count)
        arg1 = self._step_count ** (-0.5)
        arg2 = self._step_count * (self.warmup_steps ** -1.5)
        dynamic_lr = (self.d_model ** (-0.5)) * min(arg1, arg2)
        # print('dynamic_lr:', dynamic_lr)
        return [dynamic_lr for group in self.optimizer.param_groups]


def train_and_test_transformer(X_train, X_test, y_train, y_test, hyperparameters, transformer=None):
    """train transformer model on training data and evaluate model on both training and test data

    :param X_train:
    :param X_test:
    :param y_train:
    :param y_test:
    :param hyperparameters:
    :param transformer:
    :return:
    """
    print('>' * 10 + 'Hyperparameters' + '>' * 10)
    for k, v in hyperparameters.items():
        print(f'{k} = {v}')
    print('<' * 10 + 'Hyperparameters' + '<' * 10)
        
    device = hyperparameters['device']

    if transformer is None or not hyperparameters['transfer']:
        transformer = TransformerEncoderModel(ntoken=4, embedding_size=hyperparameters['embedding_size'],
                                              hidden_size=hyperparameters['hidden_size'],
                                              hidden_size_fully=hyperparameters['hidden_size_fully'],
                                              output_size=hyperparameters['output_size'],
                                              nhead=hyperparameters['nhead'],
                                              num_encoder_layers=hyperparameters['num_encoder_layers'],
                                              dropout=hyperparameters['drop_out'], softmax_bool=False,
                                              other_size=hyperparameters['other_size']).to(device)

    if hyperparameters['freezing']:
        for param in transformer.parameters():
            param.requires_grad = False
        for param in transformer.fully_connected_layers.parameters():
            param.requires_grad = True

    # if torch.cuda.device_count() > 1:
    #     print("Let's use", torch.cuda.device_count(), "GPUs!")
    #     transformer = nn.DataParallel(transformer)

    criterion = torch.nn.MSELoss()
    # optimizer = optim.SGD(rnn.parameters(), lr=0.5, momentum=0.9)
    # optimizer = torch.optim.RMSprop(rnn.parameters(), lr=0.01)
    optimizer = torch.optim.Adam(transformer.parameters(), lr=hyperparameters['lr'],
                                 weight_decay=hyperparameters['weight_decay'])
    # scheduler = CustomSchedule(optimizer, hyperparameters['embedding_size'], warm_steps=10)
    # scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.9)
    # scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lambda epoch: 1/(epoch+1))
    # scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=1)  # better to not use scheduler, gamma=1 means not using

    transformer = train_transformer(transformer, criterion, optimizer, scheduler, X_train, X_test, y_train, y_test,
                                    hyperparameters['epoch_num'], hyperparameters['batch_size'], device,
                                    hyperparameters['best_epoch'])

    return transformer


def train_and_test_transformer_order3(X_train, X_test, y_train, y_test, hyperparameters, transformer=None):
    """train transformer model on training data and evaluate model on both training and test data

    :param X_train:
    :param X_test:
    :param y_train:
    :param y_test:
    :param hyperparameters:
    :param transformer:
    :return:
    """
    print('>' * 10 + 'Hyperparameters' + '>' * 10)
    for k, v in hyperparameters.items():
        print(f'{k} = {v}')
    print('<' * 10 + 'Hyperparameters' + '<' * 10)

    device = hyperparameters['device']

    if transformer is None or not hyperparameters['transfer']:
        # transformer = TransformerEncoderModelOrder3(ntokens=[4, 16, 64], embedding_size=hyperparameters['embedding_size'],
        #                                             hidden_size=hyperparameters['hidden_size'],
        #                                             hidden_size_fully=hyperparameters['hidden_size_fully'],
        #                                             output_size=hyperparameters['output_size'],
        #                                             nhead=hyperparameters['nhead'],
        #                                             num_encoder_layers=hyperparameters['num_encoder_layers'],
        #                                             dropout=hyperparameters['drop_out'], softmax_bool=False,
        #                                             other_size=hyperparameters['other_size']).to(device)
        transformer = TransformerEncoderDecoderModelOrder3(ntokens=[4, 16, 64],
                                                           embedding_size=hyperparameters['embedding_size'],
                                                           hidden_size=hyperparameters['hidden_size'],
                                                           hidden_size_fully=hyperparameters['hidden_size_fully'],
                                                           output_size=hyperparameters['output_size'],
                                                           nhead=hyperparameters['nhead'],
                                                           num_encoder_layers=hyperparameters['num_encoder_layers'],
                                                           dropout=hyperparameters['drop_out'], softmax_bool=False,
                                                           other_size=hyperparameters['other_size']).to(device)

    if hyperparameters['freezing']:
        for param in transformer.parameters():
            param.requires_grad = False
        for param in transformer.fully_connected_layers.parameters():
            param.requires_grad = True

    # if torch.cuda.device_count() > 1:
    #     print("Let's use", torch.cuda.device_count(), "GPUs!")
    #     transformer = nn.DataParallel(transformer)

    criterion = torch.nn.MSELoss()
    # optimizer = optim.SGD(rnn.parameters(), lr=0.5, momentum=0.9)
    # optimizer = torch.optim.RMSprop(rnn.parameters(), lr=0.01)
    optimizer = torch.optim.Adam(transformer.parameters(), lr=hyperparameters['lr'],
                                 weight_decay=hyperparameters['weight_decay'])
    # scheduler = CustomSchedule(optimizer, hyperparameters['embedding_size'], warm_steps=10)
    # scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.9)
    # scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=lambda epoch: 1/(epoch+1))
    # scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=1)  # better to not use scheduler, gamma=1 means not using

    transformer = train_transformer_order3(transformer, criterion, optimizer, scheduler, X_train, X_test, y_train, y_test,
                                    hyperparameters['epoch_num'], hyperparameters['batch_size'], device,
                                    hyperparameters['best_epoch'])

    return transformer


class CombinedModel(nn.Module):
    """."""

    def __init__(self, model_list, hidden_size_fully=[10], dropout=0.1, softmax_bool=False):
        super(CombinedModel, self).__init__()
        self.model_type = 'Combined'
        self.models = nn.ModuleList(model_list)
        self.output_size = 1
        self.softmax_bool = softmax_bool

        if hidden_size_fully is None:
            self.fully_connected_layers = nn.Linear(len(self.models), self.output_size)
        elif len(hidden_size_fully) == 1:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(len(self.models), hidden_size_fully[0]),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], self.output_size), )
        elif len(hidden_size_fully) == 2:
            self.fully_connected_layers = nn.Sequential(
                nn.Linear(len(self.models), hidden_size_fully[0]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[0]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[0], hidden_size_fully[1]),
                nn.ReLU(),
                # nn.LayerNorm(hidden_size_fully[1]),
                nn.Dropout(dropout),
                nn.Linear(hidden_size_fully[1], self.output_size), )

        self._reset_parameters()

    def _reset_parameters(self):
        r"""Initiate parameters in the transformer model."""

        for p in self.fully_connected_layers.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(self, input):
        outputs = []
        for model in self.models:
            outputs.append(model(input)[0])
        combined = torch.cat(outputs, 1)
        output = self.fully_connected_layers(combined)
        if self.softmax_bool:
            output = F.softmax(output, dim=-1)

        return output, None


def save_model(model, model_dir='Model_Trained_20D', model_name='pegRNA_Model.pt'):
    start = time.time()
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    with open(os.path.join(model_dir, model_name), 'wb') as f:
        # torch.save(model.state_dict(), f)
        torch.save(model, f)
    print(f'Saving time:{time.time() - start:.2f}s')


def load_model(device, model_dir='Model_Trained_20D', model_name='pegRNA_Model.pt'):
    start = time.time()
    assert os.path.exists(os.path.join(model_dir, model_name)), 'The model does not exist!'
    model = torch.load(os.path.join(model_dir, model_name),  map_location=device)

    model.eval()
    print(f'Loading time:{time.time() - start:.2f}s')

    return model


def extract_attention(transformer, data, output_dir='Model_Trained_20D', output_name='pegRNA_attention_weight.xlsx'):
    # interpret result
    data_X = data.iloc[:, :-1]
    data_y = data.iloc[:, -1]
    _, att_w = evaluate_model.evaluate_transformer(transformer, data_X, data_y, 1024)
    t0 = time.time()
    id2char = list('ACGT')
    columns = ['Target', 'PBS', 'RT'] + [f'Target_weight({i - 20})' for i in range(47)] + \
              [f'PBS_weight({i - 16})' for i in range(17)] + [f'RT_weight({i + 1})' for i in range(20)]
    result = pd.DataFrame({c: [] for c in columns})
    for i in range(len(data_X)):
        target = "".join([id2char[id - 1] for id in data_X['Target'][i] if id > 0])
        pbs = "".join([id2char[id - 1] for id in data_X['PBS'][i] if id > 0])
        rt = "".join([id2char[id - 1] for id in data_X['RT'][i] if id > 0])
        target_w = [w if w > 0 else "NaN" for w in att_w[0][i].cpu().numpy().tolist()]
        pbs_w = [w for w in att_w[1][i].cpu().numpy().tolist() if w > 0]
        rt_w = [w if w > 0 else "NaN" for w in att_w[2][i].cpu().numpy().tolist()]
        result.loc[i] = [target, pbs, rt] + target_w + ['NaN'] * (17 - len(pbs_w)) + pbs_w + rt_w

    t1 = time.time()
    print(f'Extracting attention time:{t1 - t0}')
    result.to_excel(os.path.join(output_dir, output_name), sheet_name='attention_weight')
    print(f'Outputing time:{time.time() - t1}')


def extract_attention_order3(transformer, data, device, output_dir='Model_Trained', output_name='pegRNA_attention_weight.xlsx'):
    # interpret result
    data_X = data.iloc[:, :-1]
    data_y = data.iloc[:, -1]
    start = time.time()
    _, att_w, _ = evaluate_model.evaluate_transformer_order3(transformer, data_X, data_y, 1024, device)
    t0 = time.time()
    print(f'Calculating attention time:{t0 - start:.2f}')

    id2char = list('ACGT')
    columns = ['Target', 'PBS', 'RT'] + [f'Target_weight({i - 20})' for i in range(47)] + \
              [f'PBS_weight({i - 16})' for i in range(17)] + [f'RT_weight({i + 1})' for i in range(20)]
    result = pd.DataFrame({c: [] for c in columns})
    for i in range(len(data_X)):
        target = "".join([id2char[id - 1] for id in data_X['Target'][i] if id > 0])
        pbs = "".join([id2char[id - 1] for id in data_X['PBS'][i] if id > 0])
        rt = "".join([id2char[id - 1] for id in data_X['RT'][i] if id > 0])
        target_w = [w if w > 0 else "NaN" for w in att_w[0][i].cpu().numpy().tolist()]
        pbs_w = [w for w in att_w[1][i].cpu().numpy().tolist() if w > 0]
        rt_w = [w if w > 0 else "NaN" for w in att_w[2][i].cpu().numpy().tolist()]
        result.loc[i] = [target, pbs, rt] + target_w + ['NaN'] * (17 - len(pbs_w)) + pbs_w + rt_w

    columns = ['Target', 'PBS', 'RT'] + [f'Target_weight({i - 20})' for i in range(47-1)] + \
              [f'PBS_weight({i - 16})' for i in range(17-1)] + [f'RT_weight({i + 1})' for i in range(20-1)]
    result1 = pd.DataFrame({c: [] for c in columns})
    for i in range(len(data_X)):
        target = "".join([id2char[id - 1] for id in data_X['Target'][i] if id > 0])
        pbs = "".join([id2char[id - 1] for id in data_X['PBS'][i] if id > 0])
        rt = "".join([id2char[id - 1] for id in data_X['RT'][i] if id > 0])
        target_w = [w if w > 0 else "NaN" for w in att_w[3][i].cpu().numpy().tolist()]
        pbs_w = [w for w in att_w[4][i].cpu().numpy().tolist() if w > 0]
        rt_w = [w if w > 0 else "NaN" for w in att_w[5][i].cpu().numpy().tolist()]
        result1.loc[i] = [target, pbs, rt] + target_w + ['NaN'] * (17-1 - len(pbs_w)) + pbs_w + rt_w

    columns = ['Target', 'PBS', 'RT'] + [f'Target_weight({i - 20})' for i in range(47 - 2)] + \
              [f'PBS_weight({i - 16})' for i in range(17 - 2)] + [f'RT_weight({i + 1})' for i in range(20 - 2)]
    result2 = pd.DataFrame({c: [] for c in columns})
    for i in range(len(data_X)):
        target = "".join([id2char[id - 1] for id in data_X['Target'][i] if id > 0])
        pbs = "".join([id2char[id - 1] for id in data_X['PBS'][i] if id > 0])
        rt = "".join([id2char[id - 1] for id in data_X['RT'][i] if id > 0])
        target_w = [w if w > 0 else "NaN" for w in att_w[6][i].cpu().numpy().tolist()]
        pbs_w = [w for w in att_w[7][i].cpu().numpy().tolist() if w > 0]
        rt_w = [w if w > 0 else "NaN" for w in att_w[8][i].cpu().numpy().tolist()]
        result2.loc[i] = [target, pbs, rt] + target_w + ['NaN'] * (17-2 - len(pbs_w)) + pbs_w + rt_w

    t1 = time.time()
    print(f'Extracting attention time:{t1 - t0:.2f}')
    writer = pd.ExcelWriter(os.path.join(output_dir, output_name))
    result.to_excel(writer, '1-mer')
    result1.to_excel(writer, '2-mer')
    result2.to_excel(writer, '3-mer')
    writer.close()
    # result.to_excel(os.path.join(output_dir, output_name), sheet_name='attention_weight')
    print(f'Outputing time:{time.time() - t1:.2f}')