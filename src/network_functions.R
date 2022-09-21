# Load output files of SparCC and convert into a dataframe in long format.
dataLoader <- function(corrFile, pvalFile, colFile) {
  corrs <- RcppCNPy::npyLoad(corrFile)
  pvals <- RcppCNPy::npyLoad(pvalFile)
  
  uc90.labels <- stringr::str_split(readr::read_file(colFile), '\n')[[1]]
  uc90.labels <- uc90.labels[unlist(lapply(uc90.labels, nzchar))]
  
  
  # convert to data frame
  corrs <- as.data.frame(corrs)
  pvals <- as.data.frame(pvals)
  
  # set row and col names
  colnames(corrs) <- rownames(corrs) <- uc90.labels
  colnames(pvals) <- rownames(pvals) <- uc90.labels
  
  # set upper triangle to be zero for corrs and 1 for pvals
  corrs[upper.tri(corrs, diag=T)] <- 0
  pvals[upper.tri(pvals, diag=T)] <- 1
  
  # melt and filter
  corrs.melted <- as.matrix(corrs) %>% reshape2::melt() %>% rename(corr=value) %>% filter(corr!=0)
  pvals.melted <- as.matrix(pvals) %>% reshape2::melt() %>% rename(p=value) 
  
  # join
  mat <- corrs.melted %>% 
    left_join(pvals.melted, by=c("Var1", "Var2")) %>%
    na.omit() %>%
    mutate(
      pd = p / 2, 
      abs_corr = abs(corr)
    )
  
  return(mat)
}

# Make a tidygraph object
make_net <- function(mat, classes, alpha, min_corr, graphFile, filter_on_classes=TRUE) {
  mat.filt <- mat %>%
    filter(pd < UQ(alpha), corr >= UQ(min_corr))
  
  if (filter_on_classes) {
    mat.filt <- mat.filt %>%
      left_join(classes, by=c('Var1'='name')) %>% 
      rename(Var1_class = ResFinder_class) %>% 
      left_join(classes, by=c('Var2' = 'name')) %>% 
      rename(Var2_class = 'ResFinder_class') %>%
      filter(Var1_class != Var2_class)
  }
  
  # make tidygraph object
  graph <- tbl_graph(edges=mat.filt, directed=FALSE) %>%
    activate(nodes) %>%
    left_join(classes, by="name") %>%
    mutate(
      pagerank = centrality_pagerank(),
      nodedegree = centrality_degree()
    )
  
  # simplify
  graph.simple <- igraph::simplify(graph)
  
  # calculate coordinates in backgone layout
  bb <- graphlayouts::layout_as_backbone(graph.simple, keep=.4)
  
  # add coordinates to graph object
  graph <- graph %>%
    activate(nodes) %>%
    mutate(x=bb$xy[, 1], y=bb$xy[, 2])
  
  # export graph
  exportGraph(g=graph, prefix=graphFile)
  
  return(graph)
  
}

# Save the tidygraph in two formats: json and rds
exportGraph <- function(g, prefix){
  # write in d3 json format
  graph.d3 <- d3r::d3_igraph(g)
  write(graph.d3, file = paste0(prefix, '.json'))
  
  # write r object to a file
  saveRDS(g, file=paste0(prefix, '.rds'))
}

# function for importing json files
importGraph <- function(filename) {
  graph.d3 <- jsonlite::fromJSON(filename)
  graph <- tidygraph::tbl_graph(nodes=graph.d3$nodes, edges=graph.d3$links)
  return(graph)
}

classes_overgroups <- list(
    'Aminoglycoside' = 'Aminoglycoside/Fluoroquinolone/Quinolone', 
    'Aminoglycoside/Fluoroquinolone/Quinolone' = 'Aminoglycoside/Fluoroquinolone/Quinolone',
    'Lincosamide' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin', 
    'Lincosamide/Macrolide' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin', 
    'Lincosamide/Macrolide/Streptogramin_B' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin', 
    'Lincosamide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin_A' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Lincosamide/Pleuromutilin/Streptogramin_A' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Macrolide' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Macrolide/Streptogramin_B' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Streptogramin_A' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Streptogramin_A/Streptogramin_B' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Streptogramin_B' = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    "Lincosamide/Macrolide/Streptogramin_A/Streptogramin_B/Tetracycline" = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    "Lincosamide/Streptogramin_A" = 'Lincosamide/Macrolide/Oxazolidinone/Phenicol/Pleuromutilin/Streptogramin',
    'Oxazolidinone/Phenicol' = 'Oxazolidinone/Phenicol/Tetracycline', 
    'Oxazolidinone/Phenicol/Tetracycline' = 'Oxazolidinone/Phenicol/Tetracycline'
  )
  
# retrieve resistance classes
get_resClasses <- function() {
    require(RMariaDB)
    avaDB <- dbConnect(MariaDB(), group='AvA')
    
    query <- "select gene as name, anno_value from ResFinder_anno where anno_type='Class'"
    res <- dbSendQuery(avaDB, query)
    classes.data <- dbFetch(res)
    dbClearResult(res)
    dbDisconnect(avaDB)
    unloadNamespace('RMariaDB')
    
    classes.data2 <- classes.data %>%
      group_by(name) %>%
      summarise('ResFinder_class' = paste(anno_value, collapse="/")) %>%
      mutate(
        ResFinder_class = case_when(
          ResFinder_class %in% names(classes_overgroups) ~ as.character(classes_overgroups[ResFinder_class]), 
          TRUE ~ ResFinder_class
        )
      )
    
    return(classes.data2)
}

# Retrieve metadata
getMetadata <- function() {
  require(RMariaDB)
  avaDB <- dbConnect(MariaDB(), group='AvA')
  
  query <- "select run_accession, host, IF( (collection_date <= '2020-01-01') AND (collection_date >= '1800-01-01'), YEAR(collection_date), NULL)  as year , CONCAT(IF(SPLIT_STRING(location, ' ', 2) LIKE 'S', '-', ''), SPLIT_STRING(location, ' ', 1)) as latitude, CONCAT(IF(SPLIT_STRING(location, ' ', 4) LIKE 'W', '-', ''), SPLIT_STRING(location, ' ', 3)) as longitude, raw_reads, trimmed_fragments, instrument_platform, country, continent from Meta_public inner join run_overview using(run_accession, sample_accession, project_accession) where mapped_ResFinder=1 and mapped_Silva=1"
  res <- dbSendQuery(avaDB, query)
  classes.data <- dbFetch(res)
  dbClearResult(res)
  dbDisconnect(avaDB)
  unloadNamespace('RMariaDB')
  return(res)
}

summarise.graph <- function(net){
  # number of nodes and edges
  n_nodes <- gorder(net)
  n_edges <- gsize(net)
  
  # number of clusters / components
  dg <- groups(components(net))
  n_clusters <- length(dg)
  
  # global clustering coefficient
  cc <-  transitivity(net, type="average")
  
  # density of graph
  edge_den <- edge_density(net, loops=F)
  
  # network density
  dens <- 2 * n_edges / (n_nodes*(n_nodes-1))
  
  res <- tibble(
    number_of_nodes = n_nodes, 
    number_of_edges = n_edges,
    global_clustering_coefficient = cc,
    network_density = dens,
    edge_density=edge_den,
    number_of_components = n_clusters
  )
  return(res)
}

make_palette <- function(n, palette=NA){
  
  if (is.na(palette)) {
    hues = seq(15, 375, length = n + 1)
    colors = hcl(h=hues, l=65, c=100)[1:n]
  } else {
    getPalette <- colorRampPalette(brewer.pal(8, palette))
    colors = as.vector(getPalette(n)) 
  }
  return(colors)
}

library(readr)
get_host <- function(filename) {
  filename <- basename(filename)
  if (str_starts(filename, 'run_uc90')) {
    host <-str_replace(filename, 'run_uc90_', '')
    host <- str_replace(host, '_sparr_0.6_keep0.4.json', '')
  } else {
    host <-str_replace(filename, 'sparcc_', '')
    host <- str_replace(host, '_0.6.json', '')
  }
  if (host == "sparr_0.6_keep0.4.json") {
    host = 'all'
  }
  return(host)
  
}

extract_neighborhood <- function(gene.regex, G) {
  node.idxs <- G %>% 
    activate(nodes) %>%
    as_tibble() %>%
    rownames_to_column() %>%
    filter(str_starts(string=name,pattern =  gene.regex))
  
  if( nrow(node.idxs) == 1) {
    G.sel <- G %>%
      convert(to_local_neighborhood, node=as.numeric(node.idxs$rowname), order=1, mode="all") %>%
      activate(nodes) %>%
      mutate(sel = case_when(str_starts(name, gene.regex) ~ "Yes", TRUE ~ "No"))
      # ggraph(layout="nicely") +
      # geom_edge_link0(aes(color=corr, width=corr, alpha=corr)) +
      # geom_node_point(aes(color=ResFinder_class, shape=sel), size=3)  +
      # geom_node_label(aes(label=name), size=2, repel = T) +
      # col_scale +
      # scale_edge_color_viridis("Correlation",direction=-1) +
      # scale_edge_width(range=c(.2, 1), guide="none") +
      # scale_edge_alpha(range=c(0.5, 0.9), guide="none") +
      # scale_shape_manual("Highlighted", values=list("Yes" = 16, "No" = 15))  +
      # ggtitle(str_to_title(host)) +
      # theme_graph(base_family="sans")
  } else {
    G.sel <- NULL
  }
  
  return(G.sel)
}
