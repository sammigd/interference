library(igraph)

#upload the data
cai <- read.table(here("cai_data/cai.all.csv"), sep = " ", header = TRUE)
load(here("cai_data/cai.adjacency.RData"))

#calc degree
cai_graph = graph_from_adjacency_matrix(A)
cai_degree = degree(cai_graph)
cai_degree_i1 = ifelse(cai_degree <1, 0, 1)
cai_degree_i2 = ifelse(cai_degree <2, 0, 1)
cai_degree_i3 = ifelse(cai_degree <3, 0, 1)
cai_degree_i4 = ifelse(cai_degree <4, 0, 1)
cai_degree_i5 = ifelse(cai_degree <5, 0, 1)
cai_degree_i7 = ifelse(cai_degree <7, 0, 1)
cai_degree_i10 = ifelse(cai_degree <10, 0, 1)

cai_btwn = betweenness(cai_graph)
#cai_btwn = cut(cai_btwn, c(min(cai_btwn)-1, quantile(cai_btwn, .5), max(cai_btwn)+1))

#plot(cai_graph)

#cai_graph0 = graph_from_adjacency_matrix(A[1:11, 1:11])
#plot(cai_graph0)


test = data.frame(dgr = cai_degree,  
                  btwn = cai_btwn, 
                  cai_degree_i1 = cai_degree_i1, 
                  cai_degree_i2 = cai_degree_i2,
                  cai_degree_i3 = cai_degree_i3,
                  cai_degree_i5 = cai_degree_i5,
                  cai_degree_i10 = cai_degree_i10)
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
         neigh = cluster) %>%
  group_by(neigh) %>%
  mutate(norm_degree = (dgr - min(dgr)) / (max(dgr) - min(dgr)),
         norm_btwn = (btwn - min(btwn)) / (max(btwn) - min(btwn))) 

