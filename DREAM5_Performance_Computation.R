# R file to calculate performance metric for DREAM 5 networks 

## It is assumed your working directory is where this file is

## Clear R console screen output
cat("\014")  

## Load required libraries
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

## Get gold standard networks from DREAM5 synapse project
gold.std.ids = c('InSilico.rc1' = 'syn2787240', 'SAureus.rc1' = 'syn2787242', 
                 'EColi.rc1' = 'syn2787243', 'SCerevisiae.rc1' = 'syn2787244',
                 'InSilico.rc2' = 'syn2787240', 'SAureus.rc2' = 'syn2787242', 
                 'EColi.rc2' = 'syn2787243', 'SCerevisiae.rc2' = 'syn2787244')

## Get gene id mapping from DREAM5 synapse project
gene.map.ids = c('InSilico.rc1' = 'syn2787224', 'SAureus.rc1' = 'syn2787228', 
                 'EColi.rc1' = 'syn2787232', 'SCerevisiae.rc1' = 'syn2787236',
                 'InSilico.rc2' = 'syn2787224', 'SAureus.rc2' = 'syn2787228', 
                 'EColi.rc2' = 'syn2787232', 'SCerevisiae.rc2' = 'syn2787236')

## Get rank consensus network from synapse
rank.cons.ids = c('InSilico.rc1' = 'syn7291746', 'SAureus.rc1' = 'syn7268441', 
                  'EColi.rc1' = 'syn7301787', 'SCerevisiae.rc1' = 'syn7301801',
                  'InSilico.rc2' = 'syn7499788', 'SAureus.rc2' = 'syn7499798', 
                  'EColi.rc2' = 'syn7499856', 'SCerevisiae.rc2' = 'syn7499888')

## Get bic network from synapse
bic.ids = c('InSilico.rc1' = 'syn7291744', 'SAureus.rc1' = 'syn7268438', 
            'EColi.rc1' = 'syn7301783', 'SCerevisiae.rc1' = 'syn7301799',
            'InSilico.rc2' = 'syn7805122', 'SAureus.rc2' = 'syn7805130', 
            'EColi.rc2' = 'syn7805155', 'SCerevisiae.rc2' = 'syn7805096')

# Calculate AUPR and AUROC of rank consensus networks
performance.rc = lapply(names(rank.cons.ids), function(x, gld.std.id, gene.map.id, rank.cons.id){
  # Get gold standard network from synapse 
  net.edge.gs = read.table(synGet(gld.std.id[x])@filePath, header=F) %>%
    filter(V3 == 1)
  colnames(net.edge.gs) = c('from','to','weight')
  
  net.vert.gs = read.table(synGet(gene.map.id[x])@filePath, header = F)
  colnames(net.vert.gs) = c('name','Gene.ID')
  
  g = igraph::graph_from_data_frame(net.edge.gs, vertices = net.vert.gs)
  
  gs.adj = igraph::as_adjacency_matrix(g, type = 'both') 
  rownames(gs.adj) = V(g)$Gene.ID
  colnames(gs.adj) = V(g)$Gene.ID
  
  # Get gold standard hubs
  gs.nd = rowSums(gs.adj, na.rm = T)
  gs.hubs = names(gs.nd)[gs.nd >= (median(gs.nd) + sd(gs.nd))]
  
  # Get rank consensus network from synapse
  adj = read.csv(synGet(rank.cons.id[x])@filePath, header = T, row.names = 1)
  rownames(adj) = gsub('\\.','-', rownames(adj))
  colnames(adj) = gsub('\\.','-', colnames(adj))
  adj = adj[rownames(gs.adj), colnames(gs.adj)]
  adj = adj + t(adj)
  
  # Get hubs 
  nd = rowSums(adj, na.rm = T)
  hubs = names(nd)[nd >= (median(nd) + sd(nd))]
  
  # Compute AUPR and AUROC
  gs.adj = gs.adj + t(gs.adj)
  gs.adj[gs.adj > 1] = 1
  
  tmp.adj = as.vector(as.matrix(adj))
  tmp.gs.adj = as.vector(as.matrix(gs.adj))
  
  pred = prediction(tmp.adj, tmp.gs.adj)
  
  # AUROC
  roc.raw = performance(pred, measure = 'tpr', x.measure = 'fpr')
  auc.raw = performance(pred, measure = 'auc')@y.values[[1]]
  
  # AUPR
  pr = performance(pred, measure = 'prec')@y.values[[1]]
  rc = performance(pred, measure = 'rec')@y.values[[1]]
  aupr.raw = pracma::trapz(rc[!is.na(rc) & !is.na(pr)], pr[!is.na(rc) & !is.na(pr)])
  
  ## Compute hub enrichment stats
  tmp = fisherEnrichment(hubs, gs.hubs, rownames(adj))
  
  return(data.frame(auc = auc.raw,
                    aupr = aupr.raw,
                    pval = tmp$pval,
                    odds.ratio = tmp$odds))
}, 
c(gold.std.ids, gold.std.ids), c(gene.map.ids, gene.map.ids), rank.cons.ids)
names(performance.rc) = names(rank.cons.ids)
performance.rc = rbindlist(performance.rc, idcol = 'Name') %>%
  dplyr::rename(hub.pval = pval, hub.or = odds.ratio, auroc = auc)
  
# Calculate performance metrics for BIC networks
performance.bic = lapply(names(bic.ids), function(x, gld.std.id, gene.map.id, bic.id){
  # Get gold standard network from synapse 
  net.edge.gs = read.table(synGet(gld.std.id[x])@filePath, header=F) %>%
    filter(V3 == 1)
  colnames(net.edge.gs) = c('from','to','weight')
  
  net.vert.gs = read.table(synGet(gene.map.id[x])@filePath, header = F)
  colnames(net.vert.gs) = c('name','Gene.ID')
  
  g = igraph::graph_from_data_frame(net.edge.gs, vertices = net.vert.gs)
  
  gs.adj = igraph::as_adjacency_matrix(g, type = 'both') 
  rownames(gs.adj) = V(g)$Gene.ID
  colnames(gs.adj) = V(g)$Gene.ID
  
  # Get gold standard hubs
  gs.nd = rowSums(gs.adj, na.rm = T)
  gs.hubs = names(gs.nd)[gs.nd > (median(gs.nd) + sd(gs.nd))]
  
  # Get bic network
  load(synGet(bic.id[x])@filePath)
  adj.bic = bicNetworks$network
  rownames(adj.bic) = gsub('\\.','-', rownames(adj.bic))
  colnames(adj.bic) = gsub('\\.','-', colnames(adj.bic))
  
  adj.bic = adj.bic[rownames(gs.adj), colnames(gs.adj)]
  adj.bic = as.matrix(adj.bic)
  adj.bic = adj.bic + t(adj.bic)
  
  # Get bic hubs
  bic.nd = rowSums(adj.bic, na.rm = T)
  bic.hubs = names(bic.nd)[bic.nd > (median(bic.nd) + sd(bic.nd))]
  
  tmp.bic.adj = as.vector(as.matrix(adj.bic))
  tmp.gs.adj = as.vector(as.matrix(gs.adj))
  
  ## Compute hub enrichment stats
  tmp = fisherEnrichment(bic.hubs, gs.hubs, rownames(adj.bic))
  
  return(data.frame(tp = sum(tmp.bic.adj[tmp.gs.adj == 1] == 1)/2,
                    fp = sum(tmp.bic.adj[tmp.gs.adj == 0] == 1)/2,
                    tn = sum(tmp.bic.adj[tmp.gs.adj == 0] == 0)/2,
                    fn = sum(tmp.bic.adj[tmp.gs.adj == 1] == 0)/2,
                    pval = tmp$pval,
                    odds.ratio = tmp$odds))
}, gold.std.ids, gene.map.ids, bic.ids)
names(performance.bic) = names(bic.ids)
performance.bic = rbindlist(performance.bic, idcol = 'Name') %>%
  dplyr::mutate(mcc = (tp*tn - fp*fn)/(sqrt((tp+fn)*(tn+fp)*(tp+fp)*(tn+fn))),
                odds = (tp+tn)/(fp + fn)) %>%
  dplyr::rename(hub.pval = pval, hub.or = odds.ratio)

# Calculate performance metrics for BIC networks
performance.bic.rc = lapply(names(bic.ids), function(x, gld.std.id, gene.map.id, rank.cons.id, bic.id){
  # Get gold standard network from synapse 
  net.edge.gs = read.table(synGet(gld.std.id[x])@filePath, header=F) %>%
    filter(V3 == 1)
  colnames(net.edge.gs) = c('from','to','weight')
  
  net.vert.gs = read.table(synGet(gene.map.id[x])@filePath, header = F)
  colnames(net.vert.gs) = c('name','Gene.ID')
  
  g = igraph::graph_from_data_frame(net.edge.gs, vertices = net.vert.gs)
  
  gs.adj = igraph::as_adjacency_matrix(g, type = 'both') 
  rownames(gs.adj) = V(g)$Gene.ID
  colnames(gs.adj) = V(g)$Gene.ID
  
  # Get gold standard hubs
  gs.nd = rowSums(gs.adj, na.rm = T)
  gs.hubs = names(gs.nd)[gs.nd > (median(gs.nd) + sd(gs.nd))]
  
  # Get rank consensus network from synapse
  adj = read.csv(synGet(rank.cons.id[x])@filePath, header = T, row.names = 1)
  rownames(adj) = gsub('\\.','-', rownames(adj))
  colnames(adj) = gsub('\\.','-', colnames(adj))
  adj = adj[rownames(gs.adj), colnames(gs.adj)]
  adj = adj + t(adj)
  
  # Get bic network
  load(synGet(bic.id[x])@filePath)
  adj.bic = bicNetworks$network
  rownames(adj.bic) = gsub('\\.','-', rownames(adj.bic))
  colnames(adj.bic) = gsub('\\.','-', colnames(adj.bic))
  
  adj.bic = adj.bic[rownames(gs.adj), colnames(gs.adj)]
  adj.bic = as.matrix(adj.bic)
  adj.bic[which(adj.bic == 1)] = data.matrix(adj)[which(adj.bic == 1)]
  adj.bic = adj.bic + t(adj.bic)
  
  # Get bic hubs
  bic.nd = rowSums(adj.bic, na.rm = T)
  bic.hubs = names(bic.nd)[bic.nd > (median(bic.nd) + sd(bic.nd))]
  
  tmp.bic.adj = as.vector(as.matrix(adj.bic))
  
  # Compute AUPR and AUROC
  gs.adj = gs.adj + t(gs.adj)
  gs.adj[gs.adj > 1] = 1
  tmp.gs.adj = as.vector(as.matrix(gs.adj))
  
  pred = prediction(tmp.bic.adj, tmp.gs.adj)
  
  # AUROC
  auroc = performance(pred, measure = 'auc')@y.values[[1]]
  
  # AUPR
  pr = performance(pred, measure = 'prec')@y.values[[1]]
  rc = performance(pred, measure = 'rec')@y.values[[1]]
  aupr = pracma::trapz(rc[!is.na(rc) & !is.na(pr)], pr[!is.na(rc) & !is.na(pr)])
  
  ## Compute hub enrichment stats
  tmp = fisherEnrichment(bic.hubs, gs.hubs, rownames(adj.bic))
  
  return(data.frame(auroc = auroc,
                    aupr = aupr, 
                    pval = tmp$pval,
                    odds.ratio = tmp$odds))
}, 
gold.std.ids, gene.map.ids, rank.cons.ids, bic.ids)
names(performance.bic.rc) = names(bic.ids)
performance.bic.rc = rbindlist(performance.bic.rc, idcol = 'Name') %>%
  dplyr::rename(hub.pval = pval, hub.or = odds.ratio)

# Get github commit link
thisRepo <- getRepo(repository = "th1vairam/metanetworkSynapse", 
                    ref="branch", 
                    refName='net_eval')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath = 'DREAM5_Performance_Computation.R')

# Store results in synapse
tmp = rbindlist(list(rankConsensus = performance.rc, 
                     BIC = performance.bic, 
                     BIC.weighted = performance.bic.rc), 
                use.name = T, fill = T, idcol = 'Model')
write.table(tmp, file = 'DREAM5_Performance.tsv', row.names = F, sep = '\t')
obj = File('DREAM5_Performance.tsv', name = 'DREAM5 Performance', parentId = 'syn7248617')
obj = synStore(obj, 
               used = as.character(c(gold.std.ids, gene.map.ids, rank.cons.ids, bic.ids)), 
               executed = thisFile, 
               activityName = 'Estimate performance for DREAM5 networks')