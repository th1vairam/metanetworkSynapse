# Build rank consensus network from individual networks

# Input
folder.id = 'syn7268432' # Change the folder id argument for source networks

# Load Libraries
library(synapseClient)
library(knitr)
library(githubr)

library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(metanetwork)

synapseLogin()

# Get source networks
net.files = synQuery(paste0('select name,id from file where parentId == "', folder.id,'"')) 
ind = setdiff(grep('.csv', net.files$file.name), grep('rankConsensusNetwork', net.files$file.name))
net.files = net.files[ind,]

all.networks <- lapply(net.files$file.id, function(id){
  net <- data.table::fread(synGet(id)@filePath, stringsAsFactors=FALSE, data.table=FALSE)
  rownames(net) <- net$V1
  net <- data.matrix(net[,-1])
})
gc()

# Build rankconsensus
rankConsensus <- metanetwork::rankConsensus2(all.networks)

# Get github commit links
thisRepo <- getRepo(repository = "th1vairam/metanetworkSynapse", 
                    ref="branch", 
                    refName='rank_cons')

thisFile1 <- getPermlink(repository = thisRepo,
                         repositoryPath = 'buildConsensus2.R')

# Write results to synapse
write.csv(rankConsensus, file = 'rankConsensusNetwork2.csv', quote=F, row.names = T)
obj = File('rankConsensusNetwork2.csv', name = 'rankConsensusNetwork2.csv', parentId = folder.id)
annotations(obj) = list(fileType = 'csv', resultsType = 'network', algorithm = 'rankconsensus')
obj = synStore(obj, used = as.character(net.files$file.id), executed = thisFile1, activityName = 'Build rank consensus')