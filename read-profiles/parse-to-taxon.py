import pandas as pd
from pathlib import Path


for dir in Path('out').glob('*'):
    gene_name = dir.name
    if not dir.is_dir(): continue
    entries = []
    for file in dir.glob('*.txt'):
        sample_name = file.stem
        with open(file) as input_file:
            lines = input_file.readlines()
            for line in lines:
                count, profile = line.strip().split(' ')
                entry_dict = {'sample': sample_name, 'count': int(count), 'species': profile}
                entries.append(entry_dict)
    
    pd.DataFrame(entries).to_csv(f"out/{gene_name}_taxon.csv", index=False)