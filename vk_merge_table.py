#!/usr/bin/env python3.6
import sys
import pandas as pd

if len(sys.argv) < 3:
    print('Usage:   vk_panda.merge.py file1 file2 col_header_file1 col_header_file2 out_file')
    print('         merge two tab sepetated file based on their columns.')
    print('')
    sys.exit(1)

df1 = pd.read_table(sys.argv[1], sep="\t", header=0)
df2 = pd.read_table(sys.argv[2], sep="\t", header=0)

df_merge = pd.merge(df1, df2, how='left', left_on=sys.argv[3], right_on=sys.argv[4])
df_merge.to_csv(sys.argv[5], sep="\t", index=False, encoding='utf-8-sig')
