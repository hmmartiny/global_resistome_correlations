## NOTE: This file is just to store lists that are defined manually. 

#Define the set of hosts that we want to look at
hosts <- c('all', 'air', 'canis_lupus', 'chicken', 'cow', 'freshwater', 'homo_sapiens', 'marine', 'mouse', 'pig', 'soil')
title_hosts <- c('All', 'Air', 'Dog', 'Chicken', 'Cow', 'Freshwater', 'Human', 'Marine', 'Mouse', 'Pig', 'Soil')

# Host labels to host groups
host2group <- c(
  'air metagenome' = 'air',
  'Gallus gallus' = 'chicken',
  'chicken gut metagenome' = 'chicken',
  'Bos taurus' = 'cow',
  'cow dung metagenome' = 'cow',
  'freshwater metagenome' = 'freshwater',
  'Homo sapiens' = 'homo_sapiens',
  'marine metagenome' = 'marine',
  'Mus musculus' = 'mouse',
  'mouse metagenome' = 'mouse',
  'mouse gut metagenome' = 'mouse',
  'Sus scrofa' = 'pig',
  'Sus scrofa domesticus' = 'pig',
  'pig metagenome' = 'pig',
  'pig gut metagenome' = 'pig',
  'soil metagenome' = 'soil',
  'Canis lupus familiaris' = 'canis_lupus'
)
group2host <- tibble::enframe(host2group) %>%
  select(2:1) %>%
  tibble::deframe()

# Define color scheme for resfinder classes
classes <- c('Rifampicin',
             'Phenicol',
             'Beta_lactam',
             'Glycopeptide',
             'Aminoglycoside',
             'Aminoglycoside/Fluoroquinolone/Quinolone',
             'Macrolide',
             'Lincosamide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
             'Macrolide/Tetracycline',
             'Quinolone',
             'Folate_pathway_antagonist',
             'Lincosamide/Macrolide/Streptogramin',
             'Lincosamide/Macrolide',
             'Fosfomycin',
             'Steroid_antibacterial',
             'Lincosamide',
             'Lincosamide/Streptogramin',
             'Lincosamide/Pleuromutilin/Streptogramin',
             'Polymyxin',
             'Aminoglycoside/Fluoroquinolone/Macrolide/Phenicol/Rifampicin/Tetracycline',
             'Lincosamide/Macrolide/Streptogramin/Tetracycline',
             'Macrolide/Streptogramin',
             'Nitroimidazole',
             'Oxazolidinone/Phenicol',
             'Tetracycline',
             'Oxazolidinone/Phenicol/Tetracycline',
             'Pleuromutilin',
             'Streptogramin'
            )

classes_colors <- setNames(pals::glasbey(n=length(classes)), sort(classes))

