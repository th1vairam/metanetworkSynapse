computePerformance <- function(gold.std.adj, computed.adj, threshold = 0.5){
  # Calculate enrichment of top network hubs
  gs.nd = rowSums(gold.std.adj)
  c.nd = rowSums(computed.adj)
  
  ind.gs = which(gs.nd > mean(gs.nd, na.rm=T) + sd(gs.nd, na.rm = T))
  ind.comp = which(c.nd > mean(c.nd, na.rm=T) + sd(c.nd, na.rm = T))
  enrich.hub = CovariateAnalysis::fisherEnrichment(ind.comp, ind.gs, 1:dim(gold.std.adj)[1])
  
  # Calculate enrichment of top betweeness hubs
  gs.g = igraph::graph_from_adjacency_matrix(gold.std.adj, mode = 'undirected', weighted = T)
  c.g = igraph::graph_from_adjacency_matrix(computed.adj, mode = 'undirected', weighted = T)
  
  gs.betweness = igraph::betweenness(gs.g)
  c.betweness = igraph::betweenness(c.g)
  
  ind.gs = which(gs.betweness > (mean(gs.betweness, na.rm=T) + sd(gs.betweness, na.rm = T)))
  ind.comp = which(c.betweness > mean(c.betweness, na.rm=T) + sd(c.betweness, na.rm = T))
  enrich.bet = CovariateAnalysis::fisherEnrichment(ind.comp, ind.gs, 1:dim(gold.std.adj)[1])
  
  # Calculate performence
  gold.std.adj = as.vector(gold.std.adj)
  computed.adj = as.vector(computed.adj)
  
  pred = ROCR::prediction(computed.adj, gold.std.adj)
  
  # roc.raw = performance(pred, measure = 'tpr', x.measure = 'fpr')
  auc = ROCR::performance(pred, measure = 'auc')@y.values[[1]]
  
  pr = ROCR::performance(pred, measure = 'prec')@y.values[[1]]
  rc = ROCR::performance(pred, measure = 'rec')@y.values[[1]]
  aupr = pracma::trapz(rc[!is.na(rc) & !is.na(pr)], pr[!is.na(rc) & !is.na(pr)])
  
  F.measure = ROCR::performance(pred, measure = 'f')
  cutoff.ind = which.max(F.measure@y.values[[1]])
  # threshold = F.measure@x.values[[1]][cutoff.ind]
  F.max = F.measure@y.values[[1]][cutoff.ind]
  
  tp = sum(computed.adj[gold.std.adj == 1] >= threshold)/2
  fp = sum(computed.adj[gold.std.adj == 0] >= threshold)/2
  tn = sum(computed.adj[gold.std.adj == 0] < threshold)/2
  fn = sum(computed.adj[gold.std.adj == 1] < threshold)/2
 
   return(data.frame(auc = auc,
                    aupr = aupr,
                    F.max = F.max,
                    tp = tp,
                    fp = fp,
                    tn = tn,
                    fn = fn,
                    hub.pval = enrich.hub$pval,
                    hub.or = enrich.hub$Odds.Ratio,
                    bet.pval = enrich.bet$pval,
                    bet.or = enrich.bet$Odds.Ratio))
}