import pandas as pd

# adopted from 
# https://bitbucket.org/genomicepidemiology/globalsewage/src/master/resfinder_homology_reduce/resfinder_homology_reduce.Rmd

# load clusters
resFinderUC90 = pd.read_csv('../data/ResFinder.uc90', header=None, sep='\t')
resFinderUC90.columns = [ "record_type", "clust_number", "seq_length","H_id2target", "strand", "not_used", "not_used2", "alignment", "query", "target"]

# Record_types are: H hit, S Centroid, C cluster record
# Only need H and S
resFinderUC90 = resFinderUC90.loc[resFinderUC90['record_type'].isin(['S', 'H'])]

# Create target names
resFinderUC90.loc[resFinderUC90['record_type'] == 'S', 'target'] = resFinderUC90.loc[resFinderUC90['record_type'] == 'S', 'query']

# Load count data in long format
df = pd.read_csv('../data/resistome_data.csv', sep='\t')

# Merge homology cluster information with count data
df2 = df.merge(resFinderUC90, left_on='refSequence', right_on='query')

# Group by cluster number
df3 = df2.groupby(['run_accession', 'clust_number']).agg({'fragmentCountAln': 'sum'}).reset_index()

# Add name representative member for each cluster
df4 = df3.merge(resFinderUC90.loc[resFinderUC90['record_type'] == 'S', ['clust_number', 'target']], on='clust_number')
df4.rename(columns={'target': 'refSequence'}, inplace=True)

# save the homology cluster count data
df4[['run_accession', 'clust_number', 'refSequence', 'fragmentCountAln']].to_csv('../data/resistome_uc90_v2.csv')