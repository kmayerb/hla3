# hla3

### Description

FOR RESEARCH USE ONLY 

This repository contains Python functions for 
inferring HLA-alleles from bulk TCR beta chain data using 
a simple weight of evidence predictor.

### Notes

This method was used in:

Mayer-Blackwell K, Schattgen S, Cohen-Lavi L, Crawford JC, Souquette A, Gaevert JA, Hertz T, Thomas PG, Bradley PH, Fiore-Gartland A. TCR meta-clonotypes for biomarker discovery with tcrdist3: quantification of public, HLA-restricted TCR biomarkers of SARS-CoV-2 infection. bioRxiv (2020, https://doi.org/10.1101/2020.12.24.424260).


### Steps

Two steps are involved:

1. Tabulate presence of diagnostic TCRs in bulk samples 
2. Weigh the relative evidence of each HLA-allele per sample

### Step 1 - tabulate diagnostic TCRs in bulk samples

Tabulate the the presence of diagnostic TCRs in each repertoire. 
This is done via an exact matching (TRBV-family,CDR3) at the amino acid level
to a set of previously identified HLA-allele enriched TCRs, identified 
in the `data/HLA_associated_TCRs.tsv` file. This 
script assumes that the data has already been formated with 
TRBV genes using the IMGT nomenclature (i.e. TRBV2*01). 

The diagnostic TCRs (`tests/HLA_associated_TCRs.tsv`) were previously identified by: 
[DeWitt et al. 2018](https://elifesciences.org/articles/38358) in 
"Human T cell receptor occurrence patterns encode immune history, genetic background, and receptor specificity." *Elife* (10.7554/eLife.38358)

Tabulate exact matches against diagnostic TCRs (for details `python hla/exact.py -h`). 

In this example we tabulate these TCRs in the bulk files held in any folder.
The `resources` flag specifies the folder on all files ending with `endswith_str`. 

```
python hla/exact.py \
    --outfile bulk_files_vs_diagnostic_TCRS_templates.tsv \
    --ncpus 6 \
    --reference data/HLA_associated_TCRs.tsv \
    --resources /Volumes/T7/Emerson/ \
    --strip_str .tsv.concise.tsv.tcrdist3.tsv \
    --sep_str $"," \
    --cols_to_match v_b_gene,cdr3_b_aa \
    --cols_to_family v_b_gene \
    --col_to_count count \
    --endswith_str .tsv.concise.tsv.tcrdist3.tsv
```

### Step 2 - weigh the relative evidence of each HLA-allele per sample

Compare strength of evidence. 

Once the script in step 1 completes a new output file will be created, in this case
`bulk_files_vs_diagnostic_TCRS_templates.tsv` 
It can be loaded as a Pandas DataFrame and used as input for 
the weight of evidence of predictor as follows:

```python
import pandas as pd
from hla.predict import weight_of_evidence
df = pd.read_csv('bulk_files_vs_diagnostic_TCRS_templates.tsv', sep = '\t')
w = weight_of_evidence(hla_hits_df = df, 
    threshold = 0.1,
    locus = 'HLA-A', 
    use_detects = True)
```

#### Key Arguments
```
`threshold` : float
    The minimum weight of evidence threshold from 0-1 used to infer that a 
    sample likely comes from a person expressing the the relevant allele. 
`use_detects` : True
    This is more robust as detection of each diagnostic TCR is treated equally 
    and count information is ignored
``` 

#### Output Columns
The resulting DataFrame returned by weight_of_evidence()
has the following columns, for instance if `locus` is set to 'HLA-A'


`sample` - sample name

`threshold` - threshold used for positivity call

`method` - detection or counts

`locus` - HLA locus considered

`hla_1` - predicted HLA alelle with strongest weight of evidence if >= threshold

`hla_2` - predicted HLA alelle with second strongest weight of evidence if >= threshold

`v1` - weight of evidence value 1 associated with hla_1

`v2` - weight of evidence value 2 associated with hla_1

`p1` - predicted HLA alelle with strongest weight of evidence 

`p2` - predicted HLA alelle with strongest weight of evidence 

`HLA-A*01:01`- weight of evidence for HLA-A*01:01 relative to all evidence across all tested HLA-A alleles

`HLA-A*02:01`- weight of evidence for HLA-A*02:01 relative to all evidence across all tested HLA-A alleles

`HLA-A*03:01`- weight of evidence for HLA-A*03:01 relative to all evidence across all tested HLA-A alleles

`HLA-A*11:01`- weight of evidence for HLA-A*11:01 relative to all evidence across all tested HLA-A alleles

`HLA-A*23:01`- weight of evidence for HLA-A*23:01 relative to all evidence across all tested HLA-A alleles

`HLA-A*24:02`- weight of evidence for HLA-A*24:02 relative to all evidence across all tested HLA-A alleles

...

`HLA-A*68:02` - weight of evidence for HLA-A*68:02 relative to all evidence across all tested HLA-A alleles


### Performance

Using a fixed fixed threshold of 0.1 achieved reasonable performance. 
Based on the current diagnostic TCR set, it currently performs reasonably 
well on many of the major HLA-A and HLA-B MHC class I alleles. 


```
In [12]: performance_summary_df.query('threshold == .1 & sens > .85 & TPs > 15')
           index      sens      spec      acur threshold  TPs  TNs FPs FNs  locus  method        F1
40   HLA-A*02:01         1  0.980645  0.989071       0.1  239  304   6   0  HLA-A  detect  0.987603
41   HLA-A*01:01  0.981818  0.994792  0.990893       0.1  162  382   2   3  HLA-A  detect  0.984802
42   HLA-A*03:01   0.92126  0.995261  0.978142       0.1  117  420   2  10  HLA-A  detect   0.95122
43   HLA-A*24:02  0.910714   0.98627  0.970856       0.1  102  431   6  10  HLA-A  detect  0.927273
44   HLA-A*11:01  0.916667         1  0.992714       0.1   44  501   0   4  HLA-A  detect  0.956522
45   HLA-A*29:02         1  0.992126  0.992714       0.1   41  504   4   0  HLA-A  detect  0.953488
46   HLA-A*26:01  0.948718  0.976471  0.974499       0.1   37  498  12   2  HLA-A  detect  0.840909
47   HLA-A*32:01  0.891892  0.982422  0.976321       0.1   33  503   9   4  HLA-A  detect  0.835443
48   HLA-A*68:01         1  0.982692  0.983607       0.1   29  511   9   0  HLA-A  detect  0.865672
49   HLA-A*31:01         1  0.996169  0.996357       0.1   27  520   2   0  HLA-A  detect  0.964286
50   HLA-A*23:01      0.96  0.996183  0.994536       0.1   24  522   2   1  HLA-A  detect  0.941176
51   HLA-A*25:01         1  0.975191  0.976321       0.1   25  511  13   0  HLA-A  detect  0.793651
286  HLA-B*07:02  0.948148   0.99759  0.985455       0.1  128  414   1   7  HLA-B  detect  0.969697
288  HLA-B*44:02  0.929412         1  0.989091       0.1   79  465   0   6  HLA-B  detect  0.963415
289  HLA-B*15:01  0.964912         1  0.996364       0.1   55  493   0   2  HLA-B  detect  0.982143
290  HLA-B*35:01  0.872727         1  0.987273       0.1   48  495   0   7  HLA-B  detect  0.932039
291  HLA-B*51:01  0.932203  0.997963  0.990909       0.1   55  490   1   4  HLA-B  detect  0.956522
292  HLA-B*44:03  0.912281  0.997972  0.989091       0.1   52  492   1   5  HLA-B  detect  0.945455
293  HLA-B*18:01      0.94         1  0.994545       0.1   47  500   0   3  HLA-B  detect  0.969072
294  HLA-B*40:01  0.976744         1  0.998182       0.1   42  507   0   1  HLA-B  detect  0.988235
295  HLA-B*27:05  0.918919         1  0.994545       0.1   34  513   0   3  HLA-B  detect  0.957746
296  HLA-B*57:01  0.857143         1  0.990909       0.1   30  515   0   5  HLA-B  detect  0.923077
297  HLA-B*14:02         1  0.998081  0.998182       0.1   29  520   1   0  HLA-B  detect  0.983051
298  HLA-B*13:02         1         1         1       0.1   21  529   0   0  HLA-B  detect         1
299  HLA-B*38:01         1  0.996198  0.996364       0.1   24  524   2   0  HLA-B  detect      0.96
300  HLA-B*35:03         1   0.99811  0.998182       0.1   21  528   1   0  HLA-B  detect  0.976744
301  HLA-B*40:02         1  0.941948  0.943636       0.1   16  503  31   0  HLA-B  detect  0.507937
303  HLA-B*49:01  0.941176  0.998124  0.996364       0.1   16  532   1   1  HLA-B  detect  0.941176
621  HLA-C*07:01  0.878981  0.959391  0.936479       0.1  138  378  16  19  HLA-C  detect   0.88746
625  HLA-C*06:02  0.879518  0.963675  0.950998       0.1   73  451  17  10  HLA-C  detect  0.843931
627  HLA-C*12:03  0.885246  0.991837  0.980036       0.1   54  486   4   7  HLA-C  detect  0.907563
629  HLA-C*01:02  0.897436  0.958984  0.954628       0.1   35  491  21   4  HLA-C  detect  0.736842
630  HLA-C*08:02         1  0.992218   0.99274       0.1   37  510   4   0  HLA-C  detect  0.948718
632  HLA-C*16:01  0.916667   0.98835  0.983666       0.1   33  509   6   3  HLA-C  detect      0.88
634  HLA-C*15:02  0.857143  0.960377  0.956443       0.1   18  509  21   3  HLA-C  detect       0.6
```

The `weight_of_evidence()` function can called across a number of thresholds. 

### Simple Example 

The example shows only 4 possible alleles but actual predictions are based on full set of alleles with diagnostic TCRs.

![Simple Example](https://user-images.githubusercontent.com/46639063/126405548-1a76e898-27c4-4ff4-b82c-262f6e125955.png)

The main rationale behind using this relatively simple method is 
to accomodate the uneven number of diagnostic features per alelle and the 
differential sequence depth per sample. 

The simple approach without feature weights also was selected to prevent a classifier from 
over-fitting the relatively small amount of HLA-genotype training data. 

As additional data become available, more sophisticated ensemble ML approaches will be incorporated. 

#### Comments

##### Performance 
![fig3_weight_of_evidence](https://user-images.githubusercontent.com/46639063/126401775-66fc0a39-9079-4947-b231-acdc70f58d2a.png)

##### Input Data By Allele
![fig1_training_data](https://user-images.githubusercontent.com/46639063/126399979-8af36f19-b867-459b-8509-c7f140d4c58f.png)

##### Input Data By Allele and Probability of Generation
![fig2_training_data](https://user-images.githubusercontent.com/46639063/126399980-03cc29fb-4389-41a1-8a7a-80f9c053c386.png)


#### Citations

This method was used in:

Mayer-Blackwell K, Schattgen S, Cohen-Lavi L, Crawford JC, Souquette A, Gaevert JA, Hertz T, Thomas PG, Bradley PH, Fiore-Gartland A. TCR meta-clonotypes for biomarker discovery with tcrdist3: quantification of public, HLA-restricted TCR biomarkers of SARS-CoV-2 infection. bioRxiv (2020, https://doi.org/10.1101/2020.12.24.424260).

It uses data and insight from:

Emerson, R. O., DeWitt, W. S., Vignali, M., Gravley, J., Hu, J. K., Osborne, E. J., ... & Robins, H. S. (2017). Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nature genetics, 49(5), 659-665.

DeWitt III, W. S., Smith, A., Schoch, G., Hansen, J. A., Matsen IV, F. A., & Bradley, P. (2018). Human T cell receptor occurrence patterns encode immune history, genetic background, and receptor specificity. Elife, 7, e38358.


#### Data

Data from Emerson et al. 2017 is publicly available [here](https://clients.adaptivebiotech.com/pub/emerson-2017-natgen).


