"""
FOR RESEARCH USE ONLY 

# Drafted Jan 6, 2021
# Updated July 20, 2021
# Seattle, WA
# Copyright (c) 2021 Koshlan Mayer-Blackwell

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

"""
This script provides the input to predict.py

This command line program efficiently performs exact 
match analysis enabling rapid tabulation of TCR clonotypes in 
many bulk samples. While the reference file 
must have a single column of strings to be matched, 
greater flexibility is provided for inputs. 

It is generally assumed that users will have 
a set of bulk TCR repertoire files with 
columns, designating TRBV, TRBJ gene usages 
and the amino acid sequence of the CDR3. 

The basic usage assumes columns: 'v_b_gene' and 'cdr3_b_aa'.

Therefore the flag  --cols_to_match v_b_gene,cdr3_b_aa 

The default reference are HLA-restricted TCRs from Emerson et al. 
later identified by DeWitt et al.. They take the form
V06,CASSPGPDRYEQYF. Thus, the script must know the 
column to convert from gene to family. 

Therefore the flag --cols_to_family v_b_gene 

Finally, this script could deal with reference files 
that represent TCRs in various forms 

1. The Important Dewitt Table Form:

V06,CASSPGPDRYEQYF

Note that our input data has columns v_b_gene (TRBV06*01) and cdr3_b_aa (CASSPGPDRYEQYF)
which must be combined with a comma for exact string matching. Moreover the TRBV06*01 name must be converted 
to V06, which is achieved with flags set to 

--cols_to_family v_b_gene 
--cols_to_match  v_b_gene,cdr3_b_aa 
--sep_str $","


Other possible string matching can be achieved, for instance:

2. IMGT bioidentity form 'TRBV06*01+TRBJ06*01+CASSPGPDRYEQYF'

if flags are set to 
--cols_to_match  v_b_gene,j_b_gene,cdr3_b_aa 
--sep_str $"+"


3. Abbreviated bioidentity form 'V06+J01+CASSPGPDRYEQYF'

if flags are set to 

--cols_to_family v_b_gene,j_b_gene 
--cols_to_match  v_b_gene,j_b_gene,cdr3_b_aa 
--sep_str $"+"

or 

4. Adaptive BioIdentity From another raw bioidentity column

IMGT bioidentity form 'TCRBV06-01+TCRBJ06-01+CASSPGPDRYEQYF'

if flags are set to 

--cols_to_match bioidentity 
--sep_str $""


Notably the script will run on all files in a <resource> folder or
only a set of files specified by a comma separated string.

## FULL EXAMPLES 

Here is a full example, for copying and pasting:

python exact.py \
    --outfile test.out.tsv \
    --ncpus 6 \
    --reference data/HLA_associated_TCRs.tsv \
    --resources /Volumes/T7/immunoSEQhsTCRBV4b/tcrdist3ready \
    --strip_str .tsv.tcrdist3.tsv \
    --sep_str $"," \
    --cols_to_match v_b_gene,cdr3_b_aa \
    --cols_to_family v_b_gene

Here we only processe two files with --filenames

python exact.py \
    --outfile test.out2.tsv \
    --ncpus 6 \
    --reference data/HLA_associated_TCRs.tsv \
    --resources /Volumes/T7/immunoSEQhsTCRBV4b/tcrdist3ready \
    --filenames Subject_100.tsv.tcrdist3.tsv,Subject_101.tsv.tcrdist3.tsv \
    --strip_str .tsv.tcrdist3.tsv \
    --sep_str $"," \
    --cols_to_match v_b_gene,cdr3_b_aa \
    --cols_to_family v_b_gene

Here we scan the Full Emerson set. This example hundreds of samples. 

python exact.py \
    --outfile test.out.emerson.tsv \
    --ncpus 10 \
    --reference data/HLA_associated_TCRs.tsv \
    --resources /Volumes/T7/Emerson \
    --strip_str .tsv.concise.tsv.tcrdist3.tsv \
    --sep_str $"," \
    --cols_to_match v_b_gene,cdr3_b_aa \
    --cols_to_family v_b_gene


## OTHER OPTIONS 

Select appropriate column in input files to tabulate

--col_to_count specifies the name of the column to use for counts (which can be 'count' or 'templates' or 'productive templates')

Subset files to analyze based on suffix

--endswith_str specifies which files in resources should be analyzed (default is all .tsv)

"""
import argparse
import parmap 
import pandas as pd
import numpy as np 
import re
import os  


def get_TRV_family(s):
    """
    Converts TRBV, or TRBJ to short family represented as V[0-9]{1,2}
    
    Parameters
    ----------
    s : str

    Returns
    -------
    short : str or None

    Examples
    --------
    >>> get_TRV_family('TRBV12*01')
    'V12'
    >>> get_TRV_family('TRBV2*01')
    'V02'
    """
    try:
        r = re.search(pattern = "T[C]?R[ABGD]([VJ])([0-9]{1,2})", string = s)
        gs = r.groups()
        if len(gs[1]) == 2:
            short = f"{gs[0]}{gs[1]}"
        else:
            short = f"{gs[0]}0{gs[1]}"
    except:
        short = None
    return short

def make_str_from_list(x,sep = ','):
    """
    Parameters
    ----------
    x : list  
        list of string
    sep : str
        separator to join strings with

    Returns 
    -------
    str

    Examples 
    --------
    >>> make_str_from_list(['V1','CASAAAF'], ',')
    'V1,CASAAAF'
    >>> make_str_from_list(['V1','CASAAAF'], '+')
    'V1+CASAAAF'
    """
    return sep.join(map(str,x))

def tcrdist3_columns_to_string(df, 
    cols = ['v_b_gene','cdr3_b_aa'], 
    sep_str = ','):
    """
    Parameters
    ----------
    df : DataFrame
    
    cols : list

    sep_str : str
    
    Returns
    -------
    pd.Series of string composed of multiple columns joined by a separator string
    """
    return df[cols].apply(make_str_from_list, axis = 1, sep = sep_str )


# Do tabulation once
def t(filename, 
      resources, 
      series,
      sep = "\t",
      sep_str = ',',
      convert_to_gene_family = True,
      col_to_count = "count",
      cols_to_match = ['v_b_gene' ,'cdr3_b_aa'],
      cols_to_family = ['v_b_gene']):
    """
    tabulate 

    Parameters
    ----------

    filename : str
        File must contains cols to match
        e.g., "HIP00110.tsv.concise.tsv.tcrdist3.tsv", 
    resources : str
        e.g. destination folder= 'tests/emerson', 
    series : list of pd.Series 

    sep : str
        
    convert_to_gene_family : bool
        if True TRBV12*01 becomes V12
    cols_to_match : list
        list of columns to compose into a match string
    cols_to_family : list
        list of columns to covert from gene to family level resolution

    Notes
    -----
    1. Optionally convert columns to their gene family representation
    2. Define match column
    3. group all clones that belong to the match pattern
    4. convert all counts to a dictionary 
    5. loop through search sequences and get if in dictionary
    """
    full_path = os.path.join(resources, filename)
    df = pd.read_csv(full_path, sep = sep)
    if convert_to_gene_family:
        for col in cols_to_family:
            df[col] = df[col].apply(lambda s : get_TRV_family(s))
    df['match'] = tcrdist3_columns_to_string(df, 
        cols = cols_to_match, 
        sep_str = sep_str )
    df = df[['match', col_to_count]]
    dfg = df.groupby(['match'])[col_to_count].sum().reset_index().\
        sort_values(col_to_count, ascending = False).reset_index(drop = True)
    l = dfg.to_dict('split')['data']
    cnt ={x[0]: x[1] for x in l}
    return [cnt.get(x,0) for x in series]


# Do tabulation in parallel 
def ts( ncpus,
        filenames,
        strip_str,
        resources,
        series,
        series_hla,
        sep,
        sep_str,
        convert_to_gene_family,
        col_to_count,
        cols_to_match,
        cols_to_family):
    """
    ts is a wrapper of the function t enabled by parmap

    Parameters 
    ----------
    filenames : list 
        list of filenames in the folder <resources> to be analyzed
    ncpus : int
        how many cpus to pass to pm_processes in parmap
    
    Returns
    -------
    df : pd.DataFrame

    """
    cnts = parmap.map(t,filenames, 
        series =series, 
        resources = resources,
        sep = sep,
        col_to_count = col_to_count,
        cols_to_match = cols_to_match,
        cols_to_family = cols_to_family,
        convert_to_gene_family = convert_to_gene_family,
        pm_processes = ncpus, 
        pm_pbar = True)

    d = dict()
    fs = [f.strip(strip_str) for f in filenames]
    for k,v in zip(fs, cnts):
        d[k] = v
    df1 = pd.DataFrame({"match":series, "hla_allele": series_hla})
    df2 = pd.DataFrame(d, columns = fs)
    df = pd.concat([df1,df2], axis = 1)
    return(df)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--ncpus', 
        action="store",
        type = int,
        default = 2,
        required=False,
        help = "How many cpus to make available to parmap")
    parser.add_argument('--outfile', 
        action="store",
        type = str,
        required=True,
        help = "Where to write the final output")
    parser.add_argument('--reference', 
        action="store",
        type = str,
        default = 'data/HLA_associated_TCRs.tsv',
        required=True,
        help = "File containing HLA-diagnostic TCRs (data/HLA_associated_TCRs.tsv)")
    parser.add_argument('--resources', 
        action="store",
        type = str,
        default = 'tests/emerson',
        required=True,
        help = "path to all the files you wish to search")
    parser.add_argument('--strip_str', 
        action="store",
        type = str,
        default = '',
        required=False,
        help = "string to remove from input samples for a cleaner result")
    parser.add_argument('--endswith_str', 
        action="store",
        type = str,
        default = '.tsv',
        required=False,
        help = "What string must a file in resources directory endwith to be considered in analysis")
    parser.add_argument('--sep', 
        action="store",
        type = str,
        default = '\t',
        required=False,
        help = "This is the seperator for the input files called at the step pd.read_csv(sep = sep)")
    parser.add_argument('--sep_str', 
        action="store",
        type = str,
        default = ',' ,
        required=False,
        help = "This is the seperator between TRBV,CDR#, like , or + ")
    parser.add_argument('--cols_to_match', 
        action="store",
        type = str,
        default = 'v_b_gene,cdr3_b_aa' ,
        required=False,
        help = 'a comma seperated string like "v_b_gene,cdr3_b_aa" specifies the elements of input to form a matching string')
    parser.add_argument('--col_to_count', 
        action="store",
        type = str,
        default = 'count' ,
        required=False,
        help = "")
    parser.add_argument('--cols_to_family', 
        action="store",
        type = str,
        default = None,
        required=True,
        help = "Comma seperated string specifying features to convert for example from TRBV12*01 to V12 ")
    parser.add_argument('--filenames', 
        action="store",
        type = str,
        default = None ,
        required=False,
        help = "comma seperated list of files to run if a subset of files in resources")
    args = parser.parse_args()
    for arg in vars(args):
        print(f"{arg.upper()}={getattr(args, arg)}")
        
    
    sep                     =   args.sep
    sep_str                 =   args.sep_str
    cols_to_match           =   args.cols_to_match
    cols_to_family          =   args.cols_to_family
    col_to_count            =   args.col_to_count
    cols_to_family          =   args.cols_to_family
    
    if isinstance(cols_to_family, str):
        cols_to_family = cols_to_family.split(",")
        
    if cols_to_family is not None:
        convert_to_gene_family  =   True
    else: 
        convert_to_gene_family  =   True
    cols_to_match           =   args.cols_to_match
    
    if isinstance(cols_to_match, str):
        cols_to_match = cols_to_match.split(",")
    
    resources               =   args.resources
    filenames               =   args.filenames
    strip_str               =   args.strip_str
    endswith_str            =   args.endswith_str
    resources               =   args.resources
    ncpus                   =   args.ncpus
    outfile                 =   args.outfile

    # Load filenames, check that they are valid
    if filenames is not None:
        filenames = filenames.split(",")
        for f in filenames:
            assert os.path.isfile(os.path.join(resources, f)), f'File: {f} not found'
        print(f"RUNNING EXACT MATCH WITH {len(filenames)} VALID FILES")
    else:
        filenames = [f for f in os.listdir(resources) if f.endswith(endswith_str)]
        print(f"RUNNING EXACT MATCH WITH {len(filenames)} VALID FILES")

    # Load the reference file
    reference = pd.read_csv(args.reference, sep = sep)
    # pull the reference series of TCRs
    series      = reference['tcr']
    series_hla  = reference['hla_allele']
    
    # Do tabulation
    x = ts( ncpus = ncpus,
            filenames = filenames, 
            strip_str = strip_str,
            resources = resources, 
            series    = series,
            series_hla = series_hla,
            sep       = sep,
            sep_str   = sep_str,
            convert_to_gene_family = convert_to_gene_family,
            col_to_count           = col_to_count,
            cols_to_match          = cols_to_match ,
            cols_to_family         = cols_to_family)

    print(f"WRITING {outfile}")
    x.to_csv(outfile, sep = "\t", index = False)
    print(x)
