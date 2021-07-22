""""
FOR RESEARCH USE ONLY 

Drafted Jan 6, 2021
Updated July 20, 2021
Seattle, WA
Copyright (c) 2021 Koshlan Mayer-Blackwell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import os 
import re
import pandas as pd
import numpy as np
# ██╗    ██╗        ██████╗ ███████╗    ███████╗██╗   ██╗██╗██████╗ ███████╗███╗   ██╗ ██████╗███████╗
# ██║    ██║       ██╔═══██╗██╔════╝    ██╔════╝██║   ██║██║██╔══██╗██╔════╝████╗  ██║██╔════╝██╔════╝
# ██║ █╗ ██║       ██║   ██║█████╗      █████╗  ██║   ██║██║██║  ██║█████╗  ██╔██╗ ██║██║     █████╗  
# ██║███╗██║       ██║   ██║██╔══╝      ██╔══╝  ╚██╗ ██╔╝██║██║  ██║██╔══╝  ██║╚██╗██║██║     ██╔══╝  
# ╚███╔███╔╝██╗    ╚██████╔╝██║         ███████╗ ╚████╔╝ ██║██████╔╝███████╗██║ ╚████║╚██████╗███████╗
#  ╚══╝╚══╝ ╚═╝     ╚═════╝ ╚═╝         ╚══════╝  ╚═══╝  ╚═╝╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝╚══════╝                                                                                  
def weight_of_evidence(
    hla_hits_df,
    locus = "HLA-A",
    threshold = 0.1,
    use_detects = True,
    use_counts = False, 
    remove_columns = ['association_pvalue']):
    """

    Parameters
    ----------
    hla_hits_df : pd.DataFrame
        input DataFrame (columns are samples, rows are TCR features, values are counts per sample)
    locus : str
        "HLA-A",
    threshold : float 
        0.2
    use_detects : bool
        if True, use detections
    use_counts : bool
        if True, use counts versus detections
    remove_columns : list
        ['association_pvalue']
    
    Result 
    ------
    pd.DataFrame 
        columns:
    """
    # Remove columns that aren't feature, hla_allele, or sample>
    col_ind = [x for x in hla_hits_df.columns if x not in remove_columns]
    hla_hits_df = hla_hits_df[col_ind]
    # Subset columns to only alleles that start with <loci> string
    ind = hla_hits_df['hla_allele'].apply(lambda x : x.startswith(locus))
    hla_hits_df = hla_hits_df[ind].reset_index(drop = True)
    # Gather wide DataFrame to a Long Data Frame 
    if 'tcr' in hla_hits_df.columns:
        hla_hits_df = pd.melt(hla_hits_df, id_vars =['tcr','hla_allele'])
    if 'match' in hla_hits_df.columns:
        hla_hits_df = pd.melt(hla_hits_df, id_vars =['match','hla_allele'])

    # rebane column named variable back to sample
    hla_hits_df = hla_hits_df.rename(columns ={'variable':'sample'})
    # Summarize number of features (n) per hla_allele, (sum) of counts, and (detects)
    # Intuitively, we are looking at each sample and each allele and counting the number 
    # of diagnostic TCRs detected. 
    hla_hits_df_sum = hla_hits_df.groupby(['hla_allele','sample']).agg({'value' : ['count','sum',lambda x : np.count_nonzero(x)]}) 
    # Get ride of multi-level column names
    hla_hits_df_sum.columns = hla_hits_df_sum.columns.droplevel()
    # Get rid of row index, returning sample and allele to the dataframe
    hla_hits_df_sum = hla_hits_df_sum.reset_index()
    # Rename the columns
    hla_hits_df_sum.columns = ['hla_allele', 'sample', 'n', 'sum', 'detects']
    # <dadj> detects adjusted is detects divided by number of possible features
    # Intuitively, for each sample, and allele 
    # we are dividing the number of detects by the total possible HLA-diagnostic TCRS
    # For instance. If there were 500 possible HLA-A*02 and we detected 100 in the sample, 
    # than dadj would be 1/5. Meaning we found 20% of the diagnostic features for that allele
    hla_hits_df_sum['dadj'] = hla_hits_df_sum['detects'] /hla_hits_df_sum['n']
    # <cadj> counts adjusted is sum of counts divided by number of possible features
    hla_hits_df_sum['cadj'] = hla_hits_df_sum['sum'] /hla_hits_df_sum['n'] 
    # Compute total counts per sample, which will be left_joined below
    # Intuitvely, we now compute total adjusted counts per sample. 
    # Obviously deeper sequenced samples will potentially have more overal 
    # detects and we'll want to correct for this in the next step 
    counts_by_subject = hla_hits_df_sum.\
        groupby('sample')[['dadj','cadj']].sum().\
        reset_index()
    # Name columns once more
    counts_by_subject.columns = ['sample', 'total_dadj', 'total_cadj']

    # Left join <hla_hits_df_sum> to <counts_by_subject> so we can divide allele adjusted counts by total sample adjusted counts
    hla_hits_df_sum = hla_hits_df_sum.merge(counts_by_subject, how = "left", on = "sample")
    # <wd> weight of the evidence for allele X over total evidence, based on detects
    hla_hits_df_sum['wd'] = hla_hits_df_sum['dadj']/hla_hits_df_sum['total_dadj']
    # <wc> weight of the evidence for allele X over total evidence, based on counts
    hla_hits_df_sum['wc'] = hla_hits_df_sum['cadj']/hla_hits_df_sum['total_cadj']
    # Finally, repace NaN with 0
    hla_hits_df_sum = hla_hits_df_sum.replace(np.nan, 0)
    # Highly recommended that one uses detects
    assert use_detects != use_counts, "YOU CAN USE EITHER COUNTS (use_counts) OR DETECTS (use_detects), NOT BOTH"

    if use_detects:
        evidence = hla_hits_df_sum[['sample', 'hla_allele', 'wd']] 
    elif use_counts: 
        evidence = hla_hits_df_sum[['sample', 'hla_allele', 'wc']]

    evidence = evidence.pivot(index = ['sample'], columns = 'hla_allele')
    evidence.columns = evidence.columns.droplevel()
    # identify the alleles with the most evidence
    top2_alleles = [r.sort_values(ascending = False)[0:2].index for i,r in evidence.iterrows()]
    top2_alleles = pd.DataFrame(top2_alleles, columns = ["p1","p2"])
    top2_weights = [r.sort_values(ascending = False)[0:2].values for i,r in evidence.iterrows()]
    top2_weights = pd.DataFrame(top2_weights, columns =['v1','v2'])

    top2 = pd.concat([top2_alleles, top2_weights], axis =1 )
    # Now we apply a threshold. This is particularly necessary since the 2nd highest score is only real signal
    # if the sample comes from a heterozygous individual. 
    top2_thresholded = list()
    for i,r in top2.iterrows():
        if r['v1'] >= threshold:
            r['hla_1'] = r['p1']
        else: 
            r['hla_1'] = None
        if r['v2'] >= threshold:
            r['hla_2'] = r['p2']
        else: 
            r['hla_2'] = None
        r['threshold'] = threshold
        if use_detects: 
            r['method'] = 'detection'
        elif use_counts:
            r['method'] = 'counts'
        r['locus'] = locus

        top2_thresholded.append(r)
    top2_thresholded = pd.DataFrame(top2_thresholded)
    evidence = evidence.reset_index()
    top2_thresholded['sample'] = evidence['sample'].copy()
    # Select desired columns for final output dataframe.
    result = top2_thresholded[['sample','threshold','method','locus','hla_1','hla_2','v1','v2','p1','p2']].merge(evidence, how = "left", on = "sample")
    return result

if __name__ == "__main__":
    import pandas as pd
    import os
    from predict import weight_of_evidence
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', 
        action="store",
        type = float,
        default = 0.1,
        required=True,
        help = "")
    parser.add_argument('--input', 
        action="store",
        type = str,
        default = 'demo_files_vs_diagnostic_TCRS_templates.tsv',
        required=True,
        help = "")
    parser.add_argument('--locus', 
        action="store",
        type = str,
        default = 'HLA-A',
        required=True,
        help = "Select Locus HLA-A, HLA-B, HLA-C")
    parser.add_argument('--use_detects', 
        action="store",
        type = str,
        default = 1,
        required=False,
        help = "Use detects, set 1 to True")
    parser.add_argument('--use_counts', 
        action="store",
        type = str,
        default = 0,
        required=False,
        help = "Use counts, set 1 to True")
    parser.add_argument('--outfile', 
        action="store",
        type = str,
        default = 'test_demo_outfile.tsv',
        required=True,
        help = "filename or filepath to write predictions")

    args = parser.parse_args()
    for arg in vars(args):
        print(f"{arg.upper()}={getattr(args, arg)}")

    assert args.locus in ['HLA-A','HLA-B','HLA-C']
    assert isinstance(float(args.threshold), float)
    assert float(args.threshold) >= 0
    assert float(args.threshold) <= 1
    assert os.path.isfile(args.input)
    assert isinstance(args.outfile, str)
    
    df = pd.read_csv(args.input, sep = '\t')
    
    w = weight_of_evidence(hla_hits_df = df, 
        threshold = float(args.threshold), # 0.1
        locus = args.locus, # 'HLA-A', 
        use_detects =  bool(args.use_detects),
        use_counts  =  bool(args.use_counts))
    print(w)
    print(f"WRITING {args.outfile}")
    w.to_csv(args.outfile, sep = "\t", index = False)
