# emerson_for_tcrdist3
# 2020-12-09
import numpy as np
import pandas as pd
import os
from tcrdist.swap_gene_name import adaptive_to_imgt

def all_files(dest, suffix = ".tsv"):
    """get all files in a <dest> : str with some <suffix> : str """
    return [f for f in os.listdir(dest) if f.endswith(suffix)]

def reformat_for_tcrdist3_faster(f,d):
    df = pd.read_csv(os.path.join(d,f), sep = ',')
    print(f, df.templates.sum())
    total_templates = df['templates'].sum()
    df['subject'] = f
    df['sum_productive_templates_calc'] = total_templates
    df['cdr3_b_aa'] = df['amino_acid'].copy()
    df['count'] = df['templates'].copy()
    df['guess_v_b_gene'] = df['v_family'].apply(lambda x : f"{x}-01").apply(lambda x : adaptive_to_imgt['human'].get(x))
    df['guess_j_b_gene'] = df['j_family'].apply(lambda x : f"{x}-01").apply(lambda x : adaptive_to_imgt['human'].get(x))
    df['v_b_gene_prelim'] = df['v_gene'].apply(lambda x : adaptive_to_imgt['human'].get(x))
    df['j_b_gene_prelim'] = df['j_gene'].apply(lambda x : adaptive_to_imgt['human'].get(x))
    df['v_b_gene'] = df['v_b_gene_prelim'].copy()
    df['j_b_gene'] = df['j_b_gene_prelim'].copy()
    df.loc[df['v_b_gene_prelim'].isna(), 'v_b_gene'] = df.loc[df['v_b_gene_prelim'].isna(), 'guess_v_b_gene']
    df.loc[df['j_b_gene_prelim'].isna(), 'j_b_gene'] = df.loc[df['j_b_gene_prelim'].isna(), 'guess_j_b_gene']
    return(df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene','subject', 'count', 'productive_frequency','sum_productive_templates_calc']])

# Data from Emerson et al. 2017 is [publicly available](https://clients.adaptivebiotech.com/pub/emerson-2017-natgen)
where_are_the_files  = ''
dout = '/Volumes/T7/Emerson'
fs = all_files(dest = where_are_the_files , suffix = "concise.tsv")
for i,file in enumerate(fs):
    print(i,file)
    if not os.path.isfile(os.path.join(dout, f"{file}.tcrdist3.tsv")):
        print(i,file)
        df_tcrdist3 = reformat_for_tcrdist3_faster(f= file, d= where_are_the_files)
        df_tcrdist3.to_csv(os.path.join(dout, f"{file}.tcrdist3.tsv") , sep = "\t", index = False)

