# df = pd.read_table(f'User_Sequence/Database/output_of_search.txt', header=0)
from pegsub.models import Final_install_info, Final_correct_info
import pandas as pd

# df_correct = pd.read_table(f'User_Sequence/Database/final_pegRNA.correct.top10.txt')
# for index, row in df_correct.iterrows():
#     c = Final_correct_info()
#     c.AlleleID = row['AlleleID']
#     c.Type = row['Type']
#     c.Chromosome = row['Chromosome']
#     c.Start = row['Start']
#     c.Stop = row['Stop']
#     c.ReferenceAllele = row['ReferenceAllele']
#     c.AlternateAllele = row['AlternateAllele']
#     c.ReferenceSequence = row['ReferenceSequence(101bp)']
#     c.PegRNAID = row['PegRNAID']
#     c.Target = row['Target(47bp)']
#     c.Spacer = row['Spacer']
#     c.PBS = row['PBS']
#     c.RT = row['RT']
#     c.NickingCoordinate = row['NickingCoordinate']
#     c.Efficiency = row['Efficiency']
#     c.PAM = row['PAM']
#     c.save()
#
df_install = pd.read_table(f'User_Sequence/Database/final_pegRNA.install.top10.txt',low_memory=False)
for index, row in df_install.iterrows():
    c = Final_install_info()
    c.AlleleID = row['AlleleID']
    c.Type = row['Type']
    c.Chromosome = row['Chromosome']
    c.Start = row['Start']
    c.Stop = row['Stop']
    c.ReferenceAllele = row['ReferenceAllele']
    c.AlternateAllele = row['AlternateAllele']
    c.ReferenceSequence = row['ReferenceSequence(101bp)']
    c.PegRNAID = row['PegRNAID']
    c.Target = row['Target(47bp)']
    c.Spacer = row['Spacer']
    c.PBS = row['PBS']
    c.RT = row['RT']
    c.NickingCoordinate = row['NickingCoordinate']
    c.Efficiency = row['Efficiency']
    c.PAM = row['PAM']
    c.save()
