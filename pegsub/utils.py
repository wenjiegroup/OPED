from pegsub.models import pegRNAInfo

import pandas as pd


df=pd.read_table(f"/")
for index,row in df.iterrows():
    b=pegRNAInfo()
    b.PegRNAID=row['PegRNAID']
    b.Target=row['Target(47bp)']
    b.Spacer=row['Spacer']
    b.PBS=row['PBS']
    b.RT=row['RT']
    b.NickingCoordinate=row['NickingCoordinate']