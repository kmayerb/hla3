
if __name__ == "__main__":
    import os 
    import re
    import pandas as pd
    import numpy as np
    from hla.predict import weight_of_evidence 
    
    # allele_specific-predictions
    #  █████╗ ██╗     ██╗     ███████╗██╗     ███████╗    ███████╗██████╗ ███████╗ ██████╗██╗███████╗██╗ ██████╗
    # ██╔══██╗██║     ██║     ██╔════╝██║     ██╔════╝    ██╔════╝██╔══██╗██╔════╝██╔════╝██║██╔════╝██║██╔════╝
    # ███████║██║     ██║     █████╗  ██║     █████╗      ███████╗██████╔╝█████╗  ██║     ██║█████╗  ██║██║     
    # ██╔══██║██║     ██║     ██╔══╝  ██║     ██╔══╝      ╚════██║██╔═══╝ ██╔══╝  ██║     ██║██╔══╝  ██║██║     
    # ██║  ██║███████╗███████╗███████╗███████╗███████╗    ███████║██║     ███████╗╚██████╗██║██║     ██║╚██████╗
    # ╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝╚══════╝╚══════╝    ╚══════╝╚═╝     ╚══════╝ ╚═════╝╚═╝╚═╝     ╚═╝ ╚═════╝                                                                                                       
    #from collections import Counter
    #[x[0] for x in Counter(truth.hla_a_1.to_list() + truth.hla_a_2.to_list()).most_common() if isinstance(x[0],str)]
    #[x[0] for x in Counter(truth.hla_b_1.to_list() + truth.hla_b_2.to_list()).most_common() if isinstance(x[0],str)]
    #[x[0] for x in Counter(truth.hla_c_1.to_list() + truth.hla_c_2.to_list()).most_common() if isinstance(x[0],str)]

    #hla_hits = pd.read_csv('data/emerson.HLA_associated_TCR_QUERY.counts.tsv', sep = "\t")
    #hla_hits.columns = [x.replace('.tsv.concise.tsv', '') for x in  hla_hits.columns]
    truth = pd.read_csv('data/emerson_665_hla_truth_strings.tsv', sep = '\t')
    hla_hits = pd.read_csv('bulk_files_vs_diagnostic_TCRS_templates.tsv', sep = "\t")
    hla_hits.columns = [x.replace('.tsv.concise', '') for x in  hla_hits.columns]
    
    thresholds= [0.01,.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, .40, .45, .5]

    use_detects = True
    alleles_a = ['HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-A*24:02', 'HLA-A*11:01', 'HLA-A*29:02', 'HLA-A*26:01', 'HLA-A*32:01', 'HLA-A*68:01', 'HLA-A*31:01', 'HLA-A*23:01', 'HLA-A*25:01', 'HLA-A*30:01', 'HLA-A*68:02', 'HLA-A*02:06', 'HLA-A*30:02', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*02:05', 'HLA-A*29:01', 'HLA-A*34:01', 'HLA-A*24:03', 'HLA-A*03:02', 'HLA-A*66:01']
    alleles_b = ['HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*44:02', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*51:01', 'HLA-B*44:03', 'HLA-B*18:01', 'HLA-B*40:01', 'HLA-B*27:05', 'HLA-B*57:01', 'HLA-B*14:02', 'HLA-B*13:02', 'HLA-B*38:01', 'HLA-B*35:03', 'HLA-B*40:02', 'HLA-B*37:01', 'HLA-B*49:01', 'HLA-B*55:01', 'HLA-B*39:01', 'HLA-B*35:02', 'HLA-B*58:01', 'HLA-B*14:01', 'HLA-B*50:01', 'HLA-B*35:08', 'HLA-B*52:01', 'HLA-B*41:01', 'HLA-B*15:18', 'HLA-B*48:01', 'HLA-B*38:02', 'HLA-B*45:01', 'HLA-B*56:01', 'HLA-B*53:01', 'HLA-B*41:02', 'HLA-B*39:06', 'HLA-B*15:17', 'HLA-B*07:05', 'HLA-B*40:06', 'HLA-B*15:03', 'HLA-B*15:07']
    alleles_c = ['HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*03:04', 'HLA-C*12:03', 'HLA-C*03:03', 'HLA-C*01:02', 'HLA-C*08:02', 'HLA-C*02:02', 'HLA-C*16:01', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-C*07:04', 'HLA-C*17:01', 'HLA-C*08:01', 'HLA-C*12:02', 'HLA-C*15:05', 'HLA-C*03:02', 'HLA-C*08:03', 'HLA-C*16:02']
    alleles_lists = [alleles_a, alleles_b, alleles_c]
    loci = ["HLA-A","HLA-B","HLA-C"]
    truth_cols = ['hla_a', 'hla_b', 'hla_c']

    perf_list = list()
    granular_predictions = list()
    for locus, truth_col, alleles in zip(loci, truth_cols, alleles_lists):
        
        print(locus, truth_col, alleles)
        
        for threshold in thresholds: 
            
            print(threshold)

            predictions = weight_of_evidence( hla_hits_df= hla_hits,
                    locus = locus,
                    threshold = threshold,
                    use_detects = use_detects,
                    use_counts = False)
            
            
            
            pred_v_truth = predictions.merge(truth, how = "left", on = "sample")
            #truth_col = "hla_a"
            
            pred_v_truth = predictions.merge(truth, how = "left", on = "sample")
            indx = (pred_v_truth[truth_col].notna())&(pred_v_truth['v1'] >0)
            pred_v_truth = pred_v_truth[indx]
            print(pred_v_truth)

            perf_dict = dict()
            for a in alleles: 
                try:
                    pos = pred_v_truth[truth_col].apply(lambda x: x.find(a) != -1)
                    
                    def fmatch(x,a):
                        if not isinstance(x,str):
                            return False
                        else:
                            return x.find(a) != -1

                    match1 = pred_v_truth['hla_1'].apply(lambda x: fmatch(x,a)) 
                    match2 = pred_v_truth['hla_2'].apply(lambda x: fmatch(x,a))
                    pred_pos = match1 | match2

                    xdf = pd.DataFrame({"pos":pos, "pred_pos":pred_pos})
                    xdf['TP'] = ((xdf['pos']==True)  & (xdf['pred_pos'] == True)).astype(int)
                    xdf['FP'] = ((xdf['pos']==False) & (xdf['pred_pos'] == True)).astype(int)
                    xdf['TN'] = ((xdf['pos']==False) & (xdf['pred_pos'] == False)).astype(int)
                    xdf['FN'] = ((xdf['pos']==True)  & (xdf['pred_pos'] == False)).astype(int)
                    xdf['allele'] = a
                    xdf['truth'] = pred_v_truth['hla_a']
                    xdf['hla_1'] = pred_v_truth['hla_1']
                    xdf['hla_2'] = pred_v_truth['hla_2']
                    xdf['v1'] = pred_v_truth['v1']
                    xdf['v2'] = pred_v_truth['v2']
                    xdf['v_allele'] = pred_v_truth[a]
                    xdf['sample'] = pred_v_truth['sample']
                    xdf['threshold'] = threshold
                    if use_detects:
                        xdf['method'] = 'detect'
                    else:
                        xdf['method'] = 'count'
                    print(xdf)
                    granular_predictions.append(xdf)
                    print(xdf)

                    # summary measures
                    perf = dict()
                    perf['sens'] = xdf['TP'].sum() / (xdf['TP'].sum() + xdf['FN'].sum())
                    perf['spec'] = xdf['TN'].sum() / (xdf['TN'].sum() + xdf['FP'].sum())
                    perf['acur'] = (xdf['TP'].sum() + xdf['TN'].sum()) / (xdf['TP'].sum() + xdf['TN'].sum() + xdf['FP'].sum() + xdf['FN'].sum())
                    perf['threshold'] = threshold
                    perf['TPs'] =  xdf['TP'].sum()
                    perf['TNs'] =  xdf['TN'].sum()
                    perf['FPs'] =  xdf['FP'].sum()
                    perf['FNs'] =  xdf['FN'].sum()
                    perf['locus'] = locus
                    
                    if use_detects:
                        perf['method'] = 'detect'
                    else:
                        perf['method'] = 'count'
                    
                    print(perf)
                    perf_dict[a] = perf
                except KeyError:
                    pass
            
            perf_list.append(pd.DataFrame(perf_dict).transpose())

    performance_summary_df = pd.concat(perf_list).reset_index()
    granular_predictions_df = pd.concat(granular_predictions).reset_index()
    
    performance_summary_df['F1'] = 2*performance_summary_df.TPs / (2*performance_summary_df.TPs + performance_summary_df.FPs + performance_summary_df.FNs)
        
    performance_summary_df.to_csv('data/2021-07-20-performance_summary_hla_predictor.tsv', sep = "\t")
    granular_predictions_df.to_csv('data/2021-07-20-granular_hla_predictor1.tsv', sep = "\t")


