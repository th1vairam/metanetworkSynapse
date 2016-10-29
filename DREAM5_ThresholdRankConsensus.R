# R file to threshold rank consensus networks 

## It is assumed your working directory is where this file is

## Clear R console screen output
cat("\014")  

## Load required libraries
library(CovariateAnalysis) # Get dev branch from devtools::install_github('th1vairam/CovariateAnalysis@dev')
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(igraph)
library(ROCR)
library(pracma)
library(Matrix)

library(synapseClient)
library(githubr)

synapseLogin()

source('./pickAdjPower.R')
source('./computePerformance.R')

## Get gold standard networks from DREAM5 synapse project
GOLD.STD.IDs = c('InSilico' = 'syn2787240', 'SAureus' = 'syn2787242', 'EColi' = 'syn2787243', 'SCerevisiae' = 'syn2787244')

## Get gene id mapping from DREAM5 synapse project
GENE.MAP.IDs = c('InSilico' = 'syn2787224', 'SAureus' = 'syn2787228', 'EColi' = 'syn2787232', 'SCerevisiae' = 'syn2787236')

## Get rank consensus network from synapse
RANK.CONS.IDs = c('InSilico' = 'syn7291746', 'SAureus' = 'syn7268441', 'EColi' = 'syn7301787', 'SCerevisiae' = 'syn7301801')

## Get bic network from synapse
BIC.IDs = c('InSilico' = 'syn7291744', 'SAureus' = 'syn7268438', 'EColi' = 'syn7301783', 'SCerevisiae' = 'syn7301799')

# Get performance of all four networks
all.net.performance = mapply(function(gld.std.id, gene.map.id, rank.cons.id, bic.id){
  # Get gold standard network from synapse 
  net.edge.gs = read.table(synGet(gld.std.id)@filePath, header = FALSE, sep = '\t') %>%
    filter(V3 == 1)
  colnames(net.edge.gs) = c('from','to','weight')
  
  net.vert.gs = read.table(synGet(gene.map.id)@filePath, header = FALSE)
  colnames(net.vert.gs) = c('name','Gene.ID')
  
  g = igraph::graph_from_data_frame(net.edge.gs, vertices = net.vert.gs)
  
  gs.adj = igraph::as_adjacency_matrix(g, type = 'both') 
  gs.adj = gs.adj + t(gs.adj)
  gs.adj[gs.adj == 2] = 1
  rownames(gs.adj) = V(g)$Gene.ID
  colnames(gs.adj) = V(g)$Gene.ID
  rownames(gs.adj) = gsub('\\.','-', rownames(gs.adj))
  colnames(gs.adj) = gsub('\\.','-', colnames(gs.adj))
  
  ## Get rank consensus network from synapse
  rc.adj = read.csv(synGet(rank.cons.id)@filePath, header = T, row.names = 1) %>%
    data.matrix()
  rc.adj = rc.adj + t(rc.adj)
  rownames(rc.adj) = gsub('\\.','-', rownames(rc.adj))
  colnames(rc.adj) = gsub('\\.','-', colnames(rc.adj))
  
  # Choose appropriate power for the adjacency matrix
  adj.power.stats = pickAdjPower(rc.adj)
  ind = which(adj.power.stats$Trunc.Exp.Adj.R.Sq > 0.8)[1]
  adjusted.adj = rc.adj ^ adj.power.stats$power[ind]
  print( adj.power.stats$power[ind])
  
  ## Get bic network
  load(synGet(bic.id)@filePath)
  bic.adj = data.matrix(bicNetworks$network)
  rownames(bic.adj) = gsub('\\.','-', rownames(bic.adj))
  colnames(bic.adj) = gsub('\\.','-', colnames(bic.adj))
  bic.adj = bic.adj[rownames(gs.adj), colnames(gs.adj)]
  bic.adj[which(bic.adj == 1)] = rc.adj[which(bic.adj == 1)]
  
  rc.adj = rc.adj[rownames(gs.adj), colnames(gs.adj)]
  adjusted.adj = adjusted.adj[rownames(gs.adj), colnames(gs.adj)]
  
  all.performance = rbindlist(list(
    rc = computePerformance(gs.adj, rc.adj),
    rc.power = computePerformance(gs.adj, adjusted.adj),
    bic = computePerformance(gs.adj, bic.adj)),
    use.names = T, fill = T, idcol = 'Network')
}, GOLD.STD.IDs, GENE.MAP.IDs, RANK.CONS.IDs, BIC.IDs, SIMPLIFY = FALSE)

all.net.performance = rbindlist(all.net.performance, idcol = 'ModelName', use.name = T, fill = T)

# Store results in synapse
write.table(all.net.performance, file = 'DREAM5_Performance.tsv', sep = '\t', quote = F, row.names = F)
obj = File('DREAM5_Performance.tsv', name = 'DREAM5 Performance', parentId = 'syn7248617')
obj = synStore(obj, 
               executed = thisFile,
               used = as.character(c(GOLD.STD.IDs, GENE.MAP.IDs, RANK.CONS.IDs, BIC.IDs)),
               activityName = 'Compute DREAM5 performance')