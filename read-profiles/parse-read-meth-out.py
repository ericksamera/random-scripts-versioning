import pandas as pd
from pathlib import Path


for dir in Path('out').glob('*'):
    if not dir.is_dir(): continue
    entries = []
    for file in dir.glob('*.txt'):
        sample_name = file.stem
        with open(file) as input_file:
            entry_dict = {'sample': sample_name,}
            lines = input_file.readlines()

            total_read_count = 0
            for line in lines:
                count, profile = line.strip().split(' ')
                total_read_count += int(count)


            for line in lines:
                count, profile = line.strip().split(' ')
                entry_dict.update({f'"{profile}"': int(count)/int(total_read_count)})
            entries.append(entry_dict)

    x = pd.DataFrame(entries)
    x = x.fillna(0)
    x.to_csv(f'out/{dir.name}.csv', index=False)
    print(x)