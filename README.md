KSA2Dv2

## 1 Introduction 

This package aims to help users identify differentially expressed transcription factors (protein level) and targets (RNA) with significantly increased/reduced levels. 
The main algorithm was adapted from KSA2Dv2.

## 2 Installation and preparation 
TFA2D can be downloaded and installed in R. Installation requires devtools as a prerequisite:

```{r}
install.packages("devtools")
```
Next, install KSA2D by:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/TFA2D")
library(TFA2D)
```

Prepare 3 input files 


* protein data: Data frame of proteome data, each row is a protein, columns are samples.  
* RNA data: Data frame RNA data, each row is a protein, columns are samples.
* Files for annotating proteins and mapping TF-targert relationships can be downloaded from https://www.grnpedia.org/trrust/

Please notice that there are specific format requirements for both protData and psiteData. 

protData.tsv (in log2 scale):

```
Gene_name   Sample1 Sample2 Sample3 ...
A0AVF1  0.8 0.4 -0.2 ...
AASS    1   0.7 -0.5 ...
AATF    -0.5    0.8 NA ...
```

rnaData.tsv (in log2 scale):

```
ID   	 Sample1 Sample2 Sample3 ...
A0AUZ9  -0.004  0.3 -0.2
A0AVT1  0.3 0.1 -0.01
A0FGR8  NA  0.03    0.1
```
Please note that protData.tsv is the default data format of tidied up FragPipe TMT-I results using the [nesviLab_scripts](https://github.com/ginnyintifa/nesviLab_scripts)


## 3 Functions 

### 3.1 TF substrate information mapping
```tf_map``` mapps the user's proteome and RNA to the known TF and target mapping information 

```{r}


tf_list_sub = tf_map(protData_filename = "/Users/ginny/My Drive/nesviLab_scripts/test_output/ncc_ratio_gene_MD_tidy_mapFiile.tsv",
  rnaData_filename = "/Users/ginny/My Drive/nonCCRCC20211020/rna_ccRCC_nCCRCC/tpm_filter_log2_zeroTosmallnumber.tsv",
  tfNetwork_filename = "/Users/ginny/My Drive/nonCCRCC20211020/cc_KSA2D_TF_2/trrust_rawdata.human.tsv",
  tf_outputName = "/Users/ginny/My Drive/TFA2D/test_output/TFA2D_tf_map.tsv")


```

Return value of this function is a list ```tf_list_sub``` and can be used in the following main function:

### 3.2 Main effect model 

```balance_sample_purity_bandwidth_tf``` performs the empirical bayes modelling using the resultant lists from the ```tf_map``` function. At the same time, if purity adjustment is wanted, an annotation file should be prepared. In the annotation file, the columna name for the samples IDs should be "caseID" and the column name for purity should be "purity". Samples should be mappable to the column names of the proteome and phospho site data files. 


```{r}

balance_sample_purity_bandwidth_tf(d1_data = tf_list_sub[[1]],
d2_data = tf_list_sub[[2]],
s1_col_name = col_cc_wgii_high,
s2_col_name = col_cc_wgii_low,
nna_cutoff = 6,    
working_dir = "/Users/ginny/My Drive/TFA2D/test_output/cc_2H/",
compare_name = "cc_wgii_tf_b2",
permute_time =200,
bandwidth_factor = 2,                                   
adjust_purity_flag = T,
annot_file = "/Users/ginny/My Drive/nonCCRCC20211020/clinical/col_annot_meta_20220526.tsv")


```


here ```col_cc_wgii_high``` and ```col_cc_wgii_low``` are the sample IDs of the two comparison groups of interest. 

Important input parameters:

* ```nna_cutoff```: the minimun number of non-missing values acceptable for each feature, should be at least 6. 
* ```permute_time```: the number of permutations needed to generate the null distribution, recommend 200.
* ```compare_name```: a tag for the output files.
* ```bandwidth_factor```: a numerical value for smoothing the empirical distribution, the larger the bandwidth, the smoother the distribution gets. Default at 1. 
* ```adjust_purity_flag```: a boolean value indicate if purity needs to be adjusted when comparing the two groups. 

Return value of this function are files deposited to the ```working_dir```. If the sample sizes of the two comparison groups are very different from each other, the comparison will be done in slices. For each slice, you will find the empirical distribution of the overall signal, the fitted null distribution and the calculated differential distribution. A text file (.tsv) detailing the results will be generated as well. 

### 3.3 Prepare for cytoscape results visualization

```tfa2d_result_cytoscape``` post-precess the results and extract TF-substrate mappings with desired cutoffs and prepare edge and node files for cytoscapte visualization. 

```{r}

tfa2d_result_cytoscape(working_dir = "/Users/ginny/My Drive/TFA2D/test_output/cc_2H/",
                       output_dir = "/Users/ginny/My Drive/TFA2D/test_output/cc_2H/test1/",
                       freq_cutoff = NULL,
                       fdr_cutoff = 0.05,
                       tf_abs_fc_cutoff = 0.2,
                       substrate_abs_fc_cutoff = 0.5,
                       tfNetwork_filename = "/Users/ginny/My Drive/nonCCRCC20211020/cc_KSA2D_TF_2/trrust_rawdata.human.tsv"
)
  
    
  
```

Return values of this function are four tsv files. 
* ```merge_all.tsv```: merged TF-substrate relationships from all the slices. 
* ```merge_sig.tsv```: the significant subset of the above file. 
* ```edge.tsv```: edges extracted from the merge table for cytoscape illustration. 
* ```node.tsv```: nodes extracted from the merge table for cytoscape illustration. 

## 4 Visualization 
Recommended settings in Cytoscape:

* size: degree
* label Font size: degree
* shape: tf_bl (diamond or ellipse)
* fill color: tf_fc
* border paint: substrate_fc
* edge transparentcy: edge_sig_bl
* edge stroke color: edge_dir 
* edge width: edge_num




