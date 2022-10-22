#!/usr/bin/env python3

"""
author: Feng Liu
"""

import pandas as pd
import torch
import time
from pegRNA_PredictingCodes import read_data


def evaluate_sl(m, X_train, X_test, y_train, y_test):
    """calculate the pearson/spearman correlation coefficient for Shallow Learning methods"""
    coefficient = {}
    y_train_pred = m.predict(X_train)
    # corr_train_pred = np.corrcoef(y_train, y_train_pred)  #pearson correlation coefficient
    y_test_pred = m.predict(X_test)
    # corr_test_pred = np.corrcoef(y_test, y_test_pred) #pearson correlation coefficient
    coefficient['pearson'] = (pd.DataFrame({'ture': y_train, 'predicted': y_train_pred}).corr(method='pearson').iloc[
                                  0, 1],
                              pd.DataFrame({'ture': y_test, 'predicted': y_test_pred}).corr(method='pearson').iloc[
                                  0, 1])
    coefficient['spearman'] = (pd.DataFrame({'ture': y_train, 'predicted': y_train_pred}).corr(method='spearman').iloc[
                                   0, 1],
                               pd.DataFrame({'ture': y_test, 'predicted': y_test_pred}).corr(method='spearman').iloc[
                                   0, 1])
    return coefficient


def evaluate_rnn(rnn, X_train, y_train, X_test, y_test):
    # test
    rnn.eval()
    with torch.no_grad():
        coefficient = {'pearson': [], 'spearman': []}
        device2 = torch.device("cpu")
        rnn.to(device2)

        input = (torch.tensor(list(X_train["Target"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_train["PBS"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_train["RT"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_train["Other"]), device=device2, dtype=torch.float32))
        hidden0 = (rnn.initHidden(len(X_train['PBS']), device2),) * 3
        outputs = rnn(input, hidden0)
        # r_train = np.corrcoef(y_train, outputs.squeeze().cpu())[0, 1]
        coefficient['pearson'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs.squeeze().cpu()})
                                      .corr(method='pearson').iloc[0, 1])
        coefficient['spearman'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs.squeeze().cpu()})
                                       .corr(method='spearman').iloc[0, 1])

        input = (torch.tensor(list(X_test["Target"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_test["PBS"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_test["RT"]), device=device2, dtype=torch.float32),
                 torch.tensor(list(X_test["Other"]), device=device2, dtype=torch.float32))
        hidden0 = (rnn.initHidden(len(X_test['PBS']), device2),) * 3
        outputs = rnn(input, hidden0)
        # r_test = np.corrcoef(y_test, outputs.squeeze().cpu())[0, 1]
        coefficient['pearson'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
                                      .corr(method='pearson').iloc[0, 1])
        coefficient['spearman'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
                                       .corr(method='spearman').iloc[0, 1])
    return coefficient


def evaluate_transformer(transformer, X_train, y_train, batch_size_test, device, verbose=True):
    # test
    transformer.eval()
    n = len(X_train)
    if verbose:
        print(f'Evaluating data: {n}')
    batch_num = n // batch_size_test + 1
    start = time.time()
    with torch.no_grad():
        coefficient = {'pearson': [], 'spearman': []}
        outputs = []
        att_weights = []
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_train.iloc[start_i:end_i, :]
            # yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            outputs = outputs + output_b

            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)

        # r_train = np.corrcoef(y_train, outputs.squeeze().cpu())[0, 1]
        coefficient['pearson'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                      .corr(method='pearson').iloc[0, 1])
        coefficient['spearman'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                       .corr(method='spearman').iloc[0, 1])
        if verbose:
            print(f'Evaluating time: {time.time() - start:.2f}s')

        # input = (torch.tensor(list(X_test["Target"]), device=device2, dtype=torch.long),
        #          torch.tensor(list(X_test["PBS"]), device=device2, dtype=torch.long),
        #          torch.tensor(list(X_test["RT"]), device=device2, dtype=torch.long))
        # outputs = transformer(input)
        # # r_test = np.corrcoef(y_test, outputs.squeeze().cpu())[0, 1]
        # coefficient['pearson'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
        #                               .corr(method='pearson').iloc[0, 1])
        # coefficient['spearman'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
        #                                .corr(method='spearman').iloc[0, 1])
    return coefficient, att_weights


def evaluate_transformer_order3(transformer, X_train, y_train, batch_size_test, device, verbose=True):
    # test
    transformer.eval()
    n = len(X_train)
    if verbose:
        print(f'Evaluating data: {n}')
    batch_num = n // batch_size_test + 1
    start = time.time()
    with torch.no_grad():
        coefficient = {'pearson': [], 'spearman': []}
        outputs = []
        att_weights = []
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_train.iloc[start_i:end_i, :]
            # yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            outputs = outputs + output_b

            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)

        # r_train = np.corrcoef(y_train, outputs.squeeze().cpu())[0, 1]
        coefficient['pearson'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                      .corr(method='pearson').iloc[0, 1])
        coefficient['spearman'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                       .corr(method='spearman').iloc[0, 1])
        if verbose:
            print(f'Evaluating time: {time.time() - start:.2f}s')

        # input = (torch.tensor(list(X_test["Target"]), device=device2, dtype=torch.long),
        #          torch.tensor(list(X_test["PBS"]), device=device2, dtype=torch.long),
        #          torch.tensor(list(X_test["RT"]), device=device2, dtype=torch.long))
        # outputs = transformer(input)
        # # r_test = np.corrcoef(y_test, outputs.squeeze().cpu())[0, 1]
        # coefficient['pearson'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
        #                               .corr(method='pearson').iloc[0, 1])
        # coefficient['spearman'].append(pd.DataFrame({'ture': y_test, 'predicted': outputs.squeeze().cpu()})
        #                                .corr(method='spearman').iloc[0, 1])
    return coefficient, att_weights, outputs


def evaluate_transformer_order_optimal(transformer, X_train, y_train, batch_size_test, device, verbose=True):
    # test
    transformer.eval()
    n = len(X_train)
    if verbose:
        print(f'Evaluating data: {n}')
    batch_num = n // batch_size_test + 1
    start = time.time()
    with torch.no_grad():
        coefficient = {'pearson': [], 'spearman': []}
        outputs = []
        att_weights = []
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_train.iloc[start_i:end_i, :]
            # yb = y_train.iloc[start_i:end_i]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            outputs = outputs + output_b

            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)

        # r_train = np.corrcoef(y_train, outputs.squeeze().cpu())[0, 1]
        coefficient['pearson'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                      .corr(method='pearson').iloc[0, 1])
        coefficient['spearman'].append(pd.DataFrame({'ture': y_train, 'predicted': outputs})
                                       .corr(method='spearman').iloc[0, 1])
        if verbose:
            print(f'Evaluating time: {time.time() - start:.2f}s')

    return coefficient, att_weights, outputs


def transformer_predictor(transformer, X_test, batch_size_test, device):
    # test
    transformer.eval()
    n = len(X_test)
    print("pre")
    batch_num = n // batch_size_test + 1
    start = time.time()
    with torch.no_grad():
        outputs = []
        att_weights = []
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_test.iloc[start_i:end_i, :]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            outputs = outputs + output_b

            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)

    print(f'Predicting time: {time.time() - start}')

    return outputs, att_weights


def transformer_predictor_order3(transformer, X_test, batch_size_test, device):
    # test
    transformer.eval()
    n = len(X_test)
    batch_num = n // batch_size_test + 1
    start = time.time()
    with torch.no_grad():
        outputs = []
        att_weights = []
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_test.iloc[start_i:end_i, :]
            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            outputs = outputs + output_b

            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)

    print(f'Predicting time: {time.time() - start}')

    return outputs, att_weights


def evaluate_trained_model(transformer, gpu_id=0):
    device = torch.device(f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu")

    data = read_data.read_data_of_for_transformer()
    index_of_train_test = 38692  # first 38692 data points were used as training data, other 4457 test data
    X_train, X_test, y_train, y_test = (data.iloc[:index_of_train_test, :-1], data.iloc[index_of_train_test:, :-1],
                                        data.iloc[:index_of_train_test, -1], data.iloc[index_of_train_test:, -1])
    data_t = read_data.read_data_for_transformer_position_and_type(flag='Type')
    index_of_train_test = 3775  # first 3775 data points were used as training data, other 403 test data
    X_train_t, X_test_t, y_train_t, y_test_t = (data_t.iloc[:index_of_train_test, :-1],
                                                data_t.iloc[index_of_train_test:, :-1],
                                                data_t.iloc[:index_of_train_test, -1],
                                                data_t.iloc[index_of_train_test:, -1])
    data_p = read_data.read_data_for_transformer_position_and_type(flag='Position')
    index_of_train_test = 1774  # first 1774 data points were used as training data, other 200 test data
    X_train_p, X_test_p, y_train_p, y_test_p = (data_p.iloc[:index_of_train_test, :-1],
                                                data_p.iloc[index_of_train_test:, :-1],
                                                data_p.iloc[:index_of_train_test, -1],
                                                data_p.iloc[index_of_train_test:, -1])
    # transformer = train_model.load_model(device, model_dir='Model_Trained_20D',
    #                                      model_name='pegRNA_Model_saved.pt')  # load model

    r1 = evaluate_transformer(transformer, X_train, y_train, 1024, device)[0]
    r2 = evaluate_transformer(transformer, X_test, y_test, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib1:\n{r}')
    r1 = evaluate_transformer(transformer, X_train_t, y_train_t, 1024, device)[0]
    r2 = evaluate_transformer(transformer, X_test_t, y_test_t, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib2 Type:\n{r}')
    r1 = evaluate_transformer(transformer, X_train_p, y_train_p, 1024, device)[0]
    r2 = evaluate_transformer(transformer, X_test_p, y_test_p, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib2 Position:\n{r}')


def evaluate_trained_model_order3(transformer, gpu_id=0):
    device = torch.device(f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu")

    data = read_data.read_data_of_for_transformer_order3()
    index_of_train_test = 38692  # first 38692 data points were used as training data, other 4457 test data
    X_train, X_test, y_train, y_test = (data.iloc[:index_of_train_test, :-1], data.iloc[index_of_train_test:, :-1],
                                        data.iloc[:index_of_train_test, -1], data.iloc[index_of_train_test:, -1])
    data_t = read_data.read_data_for_transformer_position_and_type_order3(flag='Type')
    index_of_train_test = 3775  # first 3775 data points were used as training data, other 403 test data
    X_train_t, X_test_t, y_train_t, y_test_t = (data_t.iloc[:index_of_train_test, :-1],
                                                data_t.iloc[index_of_train_test:, :-1],
                                                data_t.iloc[:index_of_train_test, -1],
                                                data_t.iloc[index_of_train_test:, -1])
    data_p = read_data.read_data_for_transformer_position_and_type_order3(flag='Position')
    index_of_train_test = 1774  # first 1774 data points were used as training data, other 200 test data
    X_train_p, X_test_p, y_train_p, y_test_p = (data_p.iloc[:index_of_train_test, :-1],
                                                data_p.iloc[index_of_train_test:, :-1],
                                                data_p.iloc[:index_of_train_test, -1],
                                                data_p.iloc[index_of_train_test:, -1])
    # transformer = train_model.load_model(device, model_dir='Model_Trained_20D',
    #                                      model_name='pegRNA_Model_saved.pt')  # load model

    r1 = evaluate_transformer_order3(transformer, X_train, y_train, 1024, device)[0]
    r2 = evaluate_transformer_order3(transformer, X_test, y_test, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib1:\n{r}')
    r1 = evaluate_transformer_order3(transformer, X_train_t, y_train_t, 1024, device)[0]
    r2 = evaluate_transformer_order3(transformer, X_test_t, y_test_t, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib2 Type:\n{r}')
    r1 = evaluate_transformer_order3(transformer, X_train_p, y_train_p, 1024, device)[0]
    r2 = evaluate_transformer_order3(transformer, X_test_p, y_test_p, 1024, device)[0]
    r = pd.concat([pd.DataFrame(r1, index=['train']), pd.DataFrame(r2, index=['test'])], 0)
    print(f'Performance on Lib2 Position:\n{r}')


def transformer_predictor_order3_file(transformer, X_test, batch_size_test, device):
    # test
    transformer.eval()
    n = len(X_test)
    batch_num = n // batch_size_test + 1
    print("batch_num", batch_num)
    start = time.time()
    with torch.no_grad():
        outputs = []
        att_weights = []
        best = {'max': 0, 'Target': 'ACGT', 'PBS': 'ACGT', 'RT': 'ACGT'}
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_test.iloc[start_i:end_i, :]

            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            print("att_weights", att_weights_b)
            # 最大值并组合
            index_max = output_b.index(max(output_b))  # 获取最大值下标
            print("最大index")
            print(index_max)
            print(max(output_b))
            best['max'] = max(output_b)  # 将最大值赋给最大值字典

            df = pd.read_table(f'Sequence/pegRNA.User.txt', header=0)  # 重新读取组合表

            dfone = df.iloc[index_max]  # 获取最大值
            best['Target'] = dfone['Target(47bp)']
            best['PBS'] = dfone['PBS']
            best['RT'] = dfone['RT']

            outputs = outputs + output_b
            # print(outputs)
            # df["efficiency"] = outputs
            # df.to_csv(f'Sequence/result.txt', sep="\t")
            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)
        print(outputs)
        df["efficiency"]=outputs
        df.to_csv(f"pegsub/static/result.txt",sep="\t")
    print(f'Predicting time: {time.time() - start}')
    return best


#update
def transformer_predictor_order3_file_update(transformer, X_test, batch_size_test, device,top_n):
    # test
    transformer.eval()
    n = len(X_test)
    batch_num = n // batch_size_test + 1
    print("batch_num", batch_num)
    start = time.time()
    with torch.no_grad():
        outputs = []
        att_weights = []
        best = {'max': 0, 'Target': 'ACGT', 'PBS': 'ACGT', 'RT': 'ACGT'}
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_test.iloc[start_i:end_i, :]

            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            print("att_weights", att_weights_b)

            outputs = outputs + output_b
            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)
        # df = pd.read_table(f'Sequence/pegRNA_update.User.txt', header=0)
        df = pd.read_table(f'Sequence/pegRNA_April_update.User.txt')
        print("ds")
        print(df.columns)
        print("ds")
        print(output_b)
        df['score'] = output_b
        print(top_n)
        top_n=int(top_n)
        best_df = df.nlargest(top_n, 'score')
        # best_df = best_df[['Target(47bp)', 'PBS', 'RT', 'score']]
        best_df=best_df[['Strand','Spacer','PAM','PBS','RT','EditToNickDistance','sgRNASpacer','sgRNAPAM','NickToNickDistance','score']]
        # best_df = best_df.rename(columns={'Target(47bp)': 'Target'})
        print(best_df)
        print(outputs)
        # df["efficiency"] = outputs
        # df.to_csv(f"pegsub/static/result_update.txt", sep="\t")
    print(f'Predicting time: {time.time() - start}')
    return best_df


def transformer_predictor_order3_file_pos(transformer, X_test, batch_size_test, device,top_n):
    # test
    transformer.eval()
    n = len(X_test)
    batch_num = n // batch_size_test + 1
    print("batch_num", batch_num)
    start = time.time()
    with torch.no_grad():
        outputs = []
        att_weights = []
        best = {'max': 0, 'Target': 'ACGT', 'PBS': 'ACGT', 'RT': 'ACGT'}
        for i in range(batch_num):
            start_i = i * batch_size_test
            end_i = start_i + batch_size_test
            xb = X_test.iloc[start_i:end_i, :]

            input = (torch.tensor(list(xb["Target"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o2"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["Target_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["PBS_o3"]), device=device, dtype=torch.long),
                     torch.tensor(list(xb["RT_o3"]), device=device, dtype=torch.long))
            # torch.tensor(list(xb["Other"]), device=device, dtype=torch.float32))
            output_b, att_weights_b = transformer(input)
            output_b = output_b.squeeze(-1).cpu().numpy().tolist()
            print("att_weights", att_weights_b)

            outputs = outputs + output_b
            # print(outputs)
            # df["efficiency"] = outputs
            # df.to_csv(f'Sequence/result.txt', sep="\t")
            if not att_weights:
                att_weights = att_weights_b
            else:
                for j in range(len(att_weights_b)):
                    att_weights[j] = torch.cat((att_weights[j], att_weights_b[j]), 0)
        print(outputs)
        # df = pd.read_table(f'Sequence/pegRNA_pos.User.txt', header=0)
        df=pd.read_table(f'Position/pegRNA_April_pos.User.txt')
        print("ds")
        print(output_b)
        df['score'] = output_b
        print(top_n)
        top_n=int(top_n)
        best_df = df.nlargest(top_n, 'score')
        # best_df = best_df[['Target(47bp)', 'PBS', 'RT', 'score']]
        # best_df = best_df.rename(columns={'Target(47bp)': 'Target'})
        best_df = best_df[['Strand', 'Spacer', 'PAM', 'PBS', 'RT', 'EditToNickDistance', 'sgRNASpacer', 'sgRNAPAM',
                           'NickToNickDistance', 'score']]
        print("best_df")
        print(best_df)
        print(outputs)
        # df["efficiency"] = outputs
        # df.to_csv(f"pegsub/static/result_pos.txt", sep="\t")

    print(f'Predicting time: {time.time() - start}')
    return best_df
