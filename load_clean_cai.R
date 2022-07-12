library(igraph)

#upload the data
cai <- read.table("cai_data/cai.all.csv", sep = " ", header = TRUE)
load("cai_data/cai.adjacency.RData")

#calc degree
cai_graph = graph_from_adjacency_matrix(A)
cai_degree = degree(cai_graph)
test = data.frame(cai_degree)
test$id = row.names(test)
cai = merge(cai, test, by = 'id', all.x = T)

df = cai %>%
  mutate(trt = ifelse(intensive == 1 & delay == 0, 1, 0),
         outcome = (takeup_survey),
         risk_averse = as.vector(scale(risk_averse)),
         disaster_prob = as.vector(scale(disaster_prob)),
         rice_inc = as.vector(scale(rice_inc)), #rice as pct of income
         disaster_loss = as.vector(scale(disaster_loss)),
         age = as.vector(scale(age)),
         ricearea_2010 = as.vector(scale(ricearea_2010)),
         cluster = as.numeric(factor((village))),
         neigh = cluster) 
