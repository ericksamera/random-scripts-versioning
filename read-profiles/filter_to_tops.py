import pandas as pd
from pathlib import Path

merge_df = pd.DataFrame({'sample': ["E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18","E19","E22","E23","E20","E21","E24","E25","E26","E27","E28","E29","E30","E31","E32","E33","E34","E35","E36","E37","E38","E39","E40","E41","E42","E43","E44","E45","E46","E47","E48","E49","E50","E51","E52","E53","E54","E55","E56","E57","E58","E59","E60","E61","E62","E63","LLN-Kune"]})

for dir in Path('out').glob('*.csv'):
    if dir.stem not in ("GNAS", "H19", "IGF2R", "KCNQ1", "MEST", "NNAT", "PEG10", "PEG3", "PLAGL1", "RTL1", "SNRPN", "SUV39H1", "TXNIP", "XIST"): continue
    if dir.stem not in ("PEG3"): continue

    x = pd.read_csv(dir)
    counts = []
    for column in x:
        if column=="sample": continue
        total = x[column].sum()
        counts.append((column, total))
    
    columns_to_keep = [key[0] for key in sorted(counts, key=lambda x: x[1], reverse=True)[:200]]
    filtered_df = x[["sample"] + columns_to_keep]
    merge_df = pd.merge(merge_df, filtered_df, on='sample', how='outer')

print(merge_df)
merge_df.to_csv("output2.csv", index=False)