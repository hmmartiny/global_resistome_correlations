# Inferring correlations for the global resistome

All count and metadata is available at https://zenodo.org/record/6919377 

## Steps in analysis
### 1. Retrieve data
using a MySQL database (see https://hmmartiny.github.io/mARG/00_Data_loading.html), use
```{bash}
mysql -e "select run_accession, refSequence, IFNULL(fragmentCountAln / (refSequence_length / 1e3), 0) as fragmentCountAln from Meta_public left join AVA_public using(run_accession, sample_accession, project_accession)" > data/resistome_data.tsv
```

### 2. Create 90% identity clusters of ARGs 
To homology reduce ResFinder sequences, use
```{bash}
usearch -cluster_fast ResFinder -id 0.9 -query_cov 0.9 -target_cov 0.9 -centroids ResFinder.UC.nr90.fa -uc ResFinder.uc90
```
Then merge the clusters together, using the script [`src/uc90_convert.py`](src/uc90_convert.py)

### 3. Run SparCC
We reimplemented SparCC to run on a GPU, but the code is simply adapted from the original to use [Cupy](https://cupy.dev/) instead of Numpy. The code is in the file ['src/run_sparcc.py'](src/run_sparcc.py)

Using cuda10.1 and cupy/10.1.0, SparCC can be run for a input data matrix like this:
```{bash}
python run_sparcc.py -v -n_iter 50 -x_iter 10 -n_perm 100 -n_piter 5 -d data/resistome_uc90_pivoted.csv
```

### 4. Construct, analyze and visualize abundance data and correlation networks
See the Rmd notebook [`notebooks/Correlations.Rmd`](notebooks/Correlations.Rmd) for the rest. 