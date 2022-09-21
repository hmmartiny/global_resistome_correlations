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
classes <- c("Aminoglycoside/Fluoroquinolone/Macrolide/Phenicol/Rifampicin/Tetracycline",
             "Aminoglycoside/Fluoroquinolone/Quinolone"                                 ,
             "Beta_lactam"                                                              ,
             "Folate_pathway_antagonist"                                                ,
             "Fosfomycin"                                                               ,
             "Glycopeptide"                                                             ,
             "Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin" ,
             "Macrolide/Tetracycline"                                                   ,
             "Nitroimidazole"                                                                  ,
             "Oxazolidinone/Phenicol/Tetracycline"                                      ,
             "Phenicol"                                                                 ,
             "Pleuromutilin"                                                            ,
             "Polymyxin"                                                                ,
             "Quinolone"                                                                ,
             "Rifampicin"                                                               ,
             "Steroid_antibacterial"                                                    ,
             "Tetracycline"    )
classes_colors <- setNames(pals::alphabet(n=length(classes)), classes)

