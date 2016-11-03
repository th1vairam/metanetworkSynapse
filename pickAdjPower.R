pickAdjPower <- function(adj, power.vect = seq(100,300,5)){
adjFit = lapply(1:length(power.vect), function(i){
  adjusted.adj = adj ^ power.vect[i]
  
  node.degree = rowSums(adjusted.adj, na.rm = T)
  prb.node.degree = density(node.degree)
  
  fit1 = lm(log10(prb.node.degree$y + 1e-9) ~ log10(prb.node.degree$x + 1e-9))
  fit2 = lm(log10(prb.node.degree$y + 1e-9) ~ log10(prb.node.degree$x + 1e-9) + I(10^log10(prb.node.degree$x + 1e-9)))

  data.frame(power = power.vect[i],
             R.Sq = summary(fit1)$r.squared,
             Adj.R.Sq = summary(fit1)$adj.r.squared,
             slope = summary(fit1)$coefficients[2,1],
             Trunc.Exp.Adj.R.Sq = summary(fit2)$adj.r.squared,
             mn.k = mean(node.degree, na.rm = TRUE),
             md.k = median(node.degree, na.rm = TRUE),
             mx.k = max(node.degree, na.rm = TRUE),
             Density = sum(node.degree)/(dim(adjusted.adj)[1] * (dim(adjusted.adj)[1] - 1)),
             Centralization = dim(adjusted.adj)[1] * (max(node.degree) - mean(node.degree))/
               ((dim(adjusted.adj)[1] - 1) * (dim(adjusted.adj)[1] - 2)),
             Heterogeneity = sqrt(dim(adjusted.adj)[1] * sum(node.degree^2)/sum(node.degree)^2 - 1))
  }) %>% data.table::rbindlist(use.names = T, fill = T)
}