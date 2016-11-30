# Build BIC network from rankconsensus network

# Input
folder.id = 'syn7268837' # Change the folder id argument for source networks
expr.id = 'syn7248651' # Change the expression id argument for relevant expression file

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

# Get rankconsensus network
net.files = synQuery(paste0('select id from file where name == "rankConsensusNetwork2.csv" and parentId == "', folder.id,'"')) 

all.networks <- lapply(net.files$file.id, function(id){
  net <- data.table::fread(synGet(id)@filePath, stringsAsFactors=FALSE, data.table=FALSE)
  rownames(net) <- net$V1
  net <- data.matrix(net[,-1])
})
gc()

# Get expression data
expr = read.table(synGet(expr.id)@filePath, header = T, row.names = 1) %>%
  data.matrix() %>%
  t

# Build BIC
bicNetworks <- metanetwork::computeBICcurve(all.networks[[1]], expr, maxEdges=1e5)
save(bicNetworks,file='bicNetworks2.rda')

# Get github commit links
thisRepo <- getRepo(repository = "th1vairam/metanetworkSynapse", 
                    ref="branch", 
                    refName='rank_cons')

thisFile1 <- getPermlink(repository = thisRepo,
                         repositoryPath = 'buildBIC2.R')

# Write results to synapse
obj = File('bicNetworks2.rda', name = 'bicNetworks2.rda', parentId = folder.id)
annotations(obj) = list(fileType = 'rda', resultsType = 'network', algorithm = 'bic')
obj = synStore(obj, used = as.character(net.files$file.id), executed = thisFile1, activityName = 'Build BIC network from rankconsensus')