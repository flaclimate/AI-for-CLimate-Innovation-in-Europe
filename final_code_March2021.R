rm(list=ls(all=TRUE)) # Removes all previously created variables
gc() # frees up memory resources

setwd("/Users/francesca/Dropbox/FEEM_CMCC/Papers/Innovation4CC")
library(bibliometrix)
library(ggplot2)
library(stats)
library(utils)
library(broom)
library(plotly)
library(stringr)
library(stringi)
library(dplyr)
library(igraph)
library(scales)
library(qgraph)
library(dplyr) #data manipulation
library(gridExtra) #to view multiple plots together
library(tidytext) #text mining
library(topicmodels) #the LDA algorithm
library(tidyr) #gather()
library(dplyr) #awesome tools
library(kableExtra) #create attractive tables
library(knitr) #simple table generator
library(ggrepel) #text and label geoms for ggplot2
library(formattable) #color tile and color bar in `kables`
library(tm) #text mining
library(wordcloud)

######################################## RESEARCH
### BIBLIOMETRIC ANALYSIS
# Uploading the baseline
# The query is "innovation and climate change"
D <- c('scopus (2).bib', 'scopus (3).bib', 'scopus (4).bib',
       'scopus (5).bib', 'scopus (6).bib', 'scopus (7).bib')
M <- convert2df(D, dbsource = "scopus", format = "bibtex") #6018 articles
newM <- duplicatedMatching(M, Field = "TI", tol = 0.95) #5556 vs initial 8355 
newM <- newM %>%
  mutate(decade =
           ifelse(newM$PY %in% 1979:1989, "AR1",
                  ifelse(newM$PY %in% 1990:1994, "AR2",
                         ifelse(newM$PY %in% 1995:2000, "AR3",
                                ifelse(newM$PY %in% 2001:2006, "AR4",
                                       ifelse(newM$PY %in% 2007:2013, "AR5",
                                              ifelse(newM$PY %in% 2014:2018, "SR1.5",
                                                     ifelse(newM$PY %in% 2019:2021, "AR6",
                                                            "NA"))))))))
set.seed(23456)
fix.contractions <- function(doc) {
  # "won't" is a special case as it does not expand to "wo not"
  doc <- gsub("won't", "will not", doc)
  doc <- gsub("can't", "can not", doc)
  doc <- gsub("n't", " not", doc)
  doc <- gsub("'ll", " will", doc)
  doc <- gsub("'re", " are", doc)
  doc <- gsub("'ve", " have", doc)
  doc <- gsub("'m", " am", doc)
  doc <- gsub("'d", " would", doc)
  doc <- gsub("decision making", "decisionmaking", doc)
  doc <- gsub("decision makers", "decisionmakers", doc)
  # 's could be 'is' or could be possessive: it has no expansion
  doc <- gsub("'s", "", doc)
  return(doc)
}

newM$AB <- sapply(newM$AB, fix.contractions)
# Delete special characters
  removeSpecialChars <- function(x) gsub("[^a-zA-Z0-9 ]", " ", x)
  newM$AB <- sapply(newM$AB, removeSpecialChars)
  #convert everything to lowercase
  newM$AB <- sapply(newM$AB, tolower) 
  str(newM$AB, nchar.max = 300)
  removenum <- function(x) gsub("[0-9]+|[[:punct:]]|\\(.*\\)", " ", x)
  newM$AB<- sapply(newM$AB, removenum)
  
  stops <- c(
    tm::stopwords("english"),
    tm::stopwords("SMART"),
    "climate", "change", "are", "will", "this", "has", "new", "paper", "can", "also", "its", "iop", "journal",
    "springer", "author", "article", "articles", "project", "research", "program") %>%
    gofastr::prep_stopwords() 
  
  library(quanteda)
  library(tidyverse)
  library(RColorBrewer)
  library(stm)
  set.seed(23456)

processed <- textProcessor(newM$AB, stem = T, customstopwords = stops, metadata = newM)
plotRemoved(processed$documents, lower.thresh = seq(1, 200, by = 100))
out <- prepDocuments(processed$documents, processed$vocab, processed$meta, lower.thresh = 10)
docs <- out$documents
vocab <- out$vocab
meta <- out$meta

#Find optimal number of topics
storage1<-searchK(docs, vocab, K = c(10,15,20, 25, 30, 35, 40), 
                 data=meta, set.seed(23456), verbose=FALSE)
plot(storage1)
storage2<-searchK(docs, vocab, K = c(26,27,28,29,30,31,32,33,34,35), 
                  data=meta, set.seed(23456), verbose=FALSE)
storage3<-searchK(docs, vocab, K = c(30,31,32,33,34, 35), 
                  data=meta, set.seed(23456), verbose=FALSE)


stmfit2<-stm(out$documents,out$vocab,K=34,
             prevalence = ~decade,
             data=out$meta,seed=24601, 
             init.type = "Spectral")
labelTopics(stmfit2)
par(bty="n",col="grey40",lwd=5)
plot.STM(stmfit2,type="summary",xlim=c(0,0.10))
par(mfrow=c(2,2))
plot.STM(stmfit2, "labels", topics = c(1,2,3), labeltype = "frex", n=40, width = 80, text.cex = 0.75)
par(bty = "L", cex = 0.80, col = "mediumvioletred", font = 2, 
    col.axis = "mediumvioletred", col.lab = "mediumvioletred")
topicQuality(stmfit2,documents=out$documents)
topic.count = 34
plot(stmfit2, type = "hist", topics = sample(1:topic.count, size = 34))

plot(stmfit2, 
     type = "labels", 
     labeltype="frex",
     n = 30, 
     topics = c(29,30,31,32), #ho the nes you like or more combined 
     text.cex = 0.80, 
     width = 120)

plot(stmfit2, 
     type = "perspectives", 
     labeltype="frex",
     n = 33, 
     topics = c(9,32), 
     text.cex = 1.2)

head(stmfit2$theta) #each row is a document and each column is the proportion of the doc contained within the column index's topic


topicNames <- labelTopics(stmfit2, n = 34)
topic <- data.frame(
  TopicNumber = 1:34,
  TopicProportions = colMeans(stmfit2$theta))

model.stm.labels <- labelTopics(stmfit2, 1:topic.count)
# Estimate the effect of other variables
prep <- estimateEffect(1:topic.count ~ decade + s(PY), stmfit2, meta = out$meta, uncertainty = "Global")
stm1effect<-estimateEffect(formula=1:34 ~ PY,
                           stmobj=stmfit2,metadata=meta)
par(mfrow=c(3,3))
for (i in seq_along(sample(1:topic.count, size = 34))) {
  plot(stm1effect, "PY", method = "continuous", topics = i, 
       main = paste0(model.stm.labels$prob[i,1:3], collapse = ", "), linecol = "slategrey", printlegend = F)
}

plot(stm1effect, "PY", method = "continuous",topics = seq(i:topic.count), 
     model = model.stm.labels, printlegend = F)
#If I want to plot only some, I specify which one
plot(stm1effect, "PY", method = "continuous",topics = seq(1:3), 
     model = model.stm.labels, printlegend = T)


#Find documents associated with topics. To do this, remove the docs removed during the textprocessor
z <- newM$TI[-processed$docs.removed]
thoughts3 <- findThoughts(stmfit2, texts=z, topics=34, n=11)$docs[[1]]
plotQuote(thoughts3, width = 80, text.cex=0.75, main = "Topic 13")
thoughts <- findThoughts(stmfit2, texts=z, topics=c(1:34), n = 11)$docs[[1]]

## Explore the topic model
td_beta <- tidy(stmfit2)
td_beta %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  mutate(topic = paste0("Topic ", topic),
         term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(term, beta, fill = as.factor(topic))) +
  geom_col(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free_y") +
  coord_flip() +
  scale_x_reordered() +
  labs(x = NULL, y = expression(beta),
       title = "Highest word probabilities for each topic",
       subtitle = "Different words are associated with different topics")

td_gamma <- tidy(stmfit2, matrix = "gamma",                    
                 document_names = newM$TI)
ggplot(td_gamma, aes(gamma, fill = as.factor(topic))) +
  geom_histogram(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ topic, ncol = 3) +
  labs(title = "Distribution of document probabilities for each topic",
       y = "Number of documents", x = expression(gamma))

library(igraph)
library(visNetwork)
mod.out.corr <- topicCorr(stmfit2, cutoff = .01)
plot(mod.out.corr)

#Extract nodes and edges
links2 <- as.matrix(mod.out.corr$posadj)
net2 <- graph_from_adjacency_matrix(links2, mode = "undirected")
net2 <- igraph::simplify(net2) 

links <- igraph::as_data_frame(net2, what="edges")
nodes <- igraph::as_data_frame(net2, what="vertices")

nodes$shape <- "dot"  
nodes$title <- paste0("Topic ", topic$TopicNumber)
nodes$label <- apply(topicNames$prob, 1, function(x) paste0(x, collapse = " \n ")) # Node label. if you want to see all the words
#Otherwise do not run 
nodes$size <- (topic$TopicProportions / max(topic$TopicProportions)) * 30
nodes$font <- "18px"
nodes$id <- as.numeric(1:34)

visNetwork(nodes, links, width="100%",  height="800px", main="Research Topics") %>% 
  visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical")) %>%
  visNodes(scaling = list(max = 10)) %>%
  visIgraphLayout(layout ="layout_with_sugiyama", smooth = T) %>%
  visInteraction(navigationButtons = T)

resNetwork <- visNetwork(nodes, links, width="100%",  height="800px", main="Research Topics") %>% 
  visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical")) %>%
  visNodes(scaling = list(max = 10)) %>%
  visEdges(arrows = 'from') %>%
  visIgraphLayout(physics = TRUE, smooth = T) %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 20) %>%
  visInteraction(navigationButtons = T) 
visSave(resNetwork, "networkResearch.html", selfcontained = TRUE)

library(corrplot)
corrplot(links2, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

links3 <- as.matrix(mod.out.corr$cor)
net2 <- graph_from_adjacency_matrix(links3, mode = "directed")
plot(net2)

M <- cor(links3)
corrplot(M, method = "color", order = "AOE", tl.col = "black", tl.cex = .7)

c4 = cluster_optimal(net2)
membership(c4)
coords = layout_with_fr(net2)

V(net2)$title <- paste0("Topic ", topic$TopicNumber)
V(net2)$member <- membership(c4)
V(net2)$size <- (topic$TopicProportions / max(topic$TopicProportions)) * 30
V(net2)$font <- "18px"
V(net2)$id <- as.numeric(1:34)
plot(net2, vertex.color=membership(c4), layout=coords)

links <- igraph::as_data_frame(net2, what="edges")
nodes <- igraph::as_data_frame(net2, what="vertices")

par(mfrow=c(1,2))
plot(net2, vertex.color=membership(c4), layout=coords)
corrplot(M, method = "color", order = "AOE", tl.col = "black", tl.cex = .7)

plot(stmfit2, 
     type = "perspectives", 
     labeltype="prob",
     n = 33, 
     topics = c(17, 24), 
     text.cex = 0.8)
## Explore similar and dissimilar topics with this plot
topicQuality(stmfit2, out$documents)


## Continuous effect: understand the time evolution
stm1effect<-estimateEffect(formula=1:33 ~ PY,
                           stmobj=stmfit2,metadata=meta) #smooth the effect of year

#plot.estimateEffect(stm1effect,         covariate= "PY", model=stmfit2,
                                        #topics=stm1effect$topics[1],
                                        #method="continuous",
                                        #xlab="Election Year",
                                        #ylab="Expected Topic Proportions",
                                        #main="Importance of ", linecol = "green",
                                        #printlegend=F)


par(mfrow=c(3,3))
for (i in seq_along(sample(1:topic.count, size = 9))) {
  plot(stm1effect, "PY", method = "continuous", topics = i, 
       main = paste0(model.stm.labels$prob[i,1:3], collapse = ", "), linecol = "slategrey", printlegend = F)
}

stm1effect<-estimateEffect(formula=1:33 ~ AU1_UN,
                           stmobj=stmfit2,metadata=meta) #smooth the effect of year
par(mfrow=c(2,2))
for (i in seq_along(sample(1:topic.count, size = 33))) {
  plot(stm1effect, "decade", method = "pointestimate", topics = i, 
       main = paste0(model.stm.labels$prob[i,1:3], collapse = ", "), linecol = "slategrey", printlegend = F)
}

#Topic prevalence in documents
plot(stmfit2, type = "hist", topics = sample(1:topic.count, size = 33))

## GRAPHICS
library(ggthemes)
library(kableExtra)
#Observe probabilities that each word is generated from each topic
td_beta <- tidy(stmfit2)
#Observe probabilities that each document is generated from each topic
td_gamma <- tidy(stmfit2, matrix = "gamma",
                 document_names = newM$TI)

# Here I am extracting the top 10 terms per topic
top_terms <- td_beta %>%
  arrange(beta) %>%
  group_by(topic) %>%
  top_n(20, beta) %>%
  arrange(-beta) %>%
  select(topic, term) %>%
  summarise(terms = list(term)) %>%
  mutate(terms = map(terms, paste, collapse = ", ")) 

gamma_terms <- td_gamma %>%
  group_by(topic) %>%
  summarise(gamma = mean(gamma)) %>%
  arrange(desc(gamma)) %>%
  left_join(top_terms, by = "topic") %>%
  mutate(topic = paste0("Topic ", topic),
         topic = reorder(topic, gamma))

gamma_terms %>%
  top_n(20, gamma) %>%
  ggplot(aes(topic, gamma, label = terms, fill = topic)) +
  geom_col(show.legend = FALSE) +
  geom_text(hjust = 0, nudge_y = 0.0005, size = 3,
            family = "IBMPlexSans") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 0.09),
                     labels = percent_format()) +
  theme_tufte(base_family = "IBMPlexSans", ticks = FALSE) +
  theme(plot.title = element_text(size = 16,
                                  family="IBMPlexSans-Bold"),
        plot.subtitle = element_text(size = 13)) +
  labs(x = NULL, y = expression(gamma),
       title = "Top 20 topics by prevalence",
       subtitle = "With the top words that contribute to each topic")

gamma_terms %>%
  select(topic, gamma, terms) %>%
  kable(digits = 3, col.names = c("Topic", "Expected topic proportion", "Top 7 terms"))

##Top documents per topic
td_gamma_top <- td_gamma %>%
  arrange(gamma) %>%
  group_by(topic) %>%
  top_n(3, gamma) %>%
  arrange(-gamma) %>%
  select(topic, document) %>%
  summarise(document = list(document)) %>% ##if you remove this line, you will have the top documents in decreasing order ordered one after the other
  mutate(document = map(document, paste, collapse = "; ")) %>%
  unnest(cols = document) 

final_docword <- td_gamma_top %>% left_join(top_terms, by = "topic")
kbl(final_docword) %>%
  kable_styling(bootstrap_options = c("condensed", "responsive"))

kbl(final_docword) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = TRUE, border_right = T) %>%
  column_spec(2, width = "30em") %>%
  save_kable(file = "table2.html", self_contained = T)

################################## PROJECTS
h2020 <- read.csv("h2020_CC.csv") #951 projects
fp7 <- read.csv("fp7_CC.csv") #918 projects
fp6 <- read.csv("fp6_CC.csv") #198 projects
projects <- read.csv("projects_CC.csv")
set.seed(23456)

projects_fp <- projects %>% count(frameworkProgramme, sort = TRUE)
projects$ecMaxContribution <- as.numeric(as.character(projects$ecMaxContribution))
projects$date <- as.Date(projects$startDate, format = "%d/%m/%Y")

ggplot(projects, aes(x=date, y=ecMaxContribution)) +
  geom_line()

fix.contractions <- function(doc) {
  # "won't" is a special case as it does not expand to "wo not"
  doc <- gsub("won't", "will not", doc)
  doc <- gsub("can't", "can not", doc)
  doc <- gsub("n't", " not", doc)
  doc <- gsub("'ll", " will", doc)
  doc <- gsub("'re", " are", doc)
  doc <- gsub("'ve", " have", doc)
  doc <- gsub("'m", " am", doc)
  doc <- gsub("'d", " would", doc)
  doc <- gsub("decision making", "decisionmaking", doc)
  doc <- gsub("decision makers", "decisionmakers", doc)
  # 's could be 'is' or could be possessive: it has no expansion
  doc <- gsub("'s", "", doc)
  return(doc)
}

projects$objective <- sapply(projects$objective, fix.contractions)
# Delete special characters
removeSpecialChars <- function(x) gsub("[^a-zA-Z0-9 ]", " ", x)
projects$objective <- sapply(projects$objective, removeSpecialChars)
#convert everything to lowercase
projects$objective <- sapply(projects$objective, tolower) 
str(projects$objective, nchar.max = 300)
removenum <- function(x) gsub("[0-9]+|[[:punct:]]|\\(.*\\)", " ", x)
projects$objective <- sapply(projects$objective, removenum)

date <- str_split_fixed(projects$startDate, "/", 3)
date <- as.data.frame(apply(date, 2, as.numeric))
date <- date[, -1]
date <- as.data.frame(date)
projects <- cbind(projects, date)

stops <- c("climate", "change", "new", "paper", "iop", "journal",
  "springer", "author", "article", "articles", "innovation", "will") %>%
  gofastr::prep_stopwords() 

set.seed(23456)
p_processed <- textProcessor(projects$objective, customstopwords = stops, stem = T, metadata = projects)
plotRemoved(p_processed$documents, lower.thresh = seq(1, 200, by = 100))
p_out <- prepDocuments(p_processed$documents, p_processed$vocab, p_processed$meta, lower.thresh = 10)
p_docs <- p_out$documents
p_vocab <- p_out$vocab
p_meta <- p_out$meta

p_storage1<-searchK(p_docs, p_vocab, K = c(10,15,20, 25, 30, 35, 40), 
                  data=p_meta, set.seed(23456), verbose=FALSE)
plot(p_storage1)
p_storage2<-searchK(p_docs, p_vocab, K = c(30,31,32,33,34,35), 
                  data=p_meta, set.seed(23456), verbose=FALSE)

p_stmfit2<-stm(p_out$documents,p_out$vocab,K=33,
             prevalence = ~frameworkProgramme,
             data=p_out$meta,seed=24601, 
             init.type = "Spectral")
labelTopics(p_stmfit2)
plot.STM(p_stmfit2,type="summary",xlim=c(0,0.10))
par(bty="n",col="grey40",lwd=5)
plot.STM(p_stmfit2,type="summary",xlim=c(0,0.10))

topicQuality(p_stmfit2,documents=p_out$documents)

p_topicNames <- labelTopics(p_stmfit2, n = 33)
p_topic <- data.frame(
  TopicNumber = 1:33,
  TopicProportions = colMeans(p_stmfit2$theta))

# get mean topic proportions per decade
topic_proportion_per_decade <- aggregate(p_stmfit2$theta, by = list(year = projects$V3), mean)

# reshape data frame
library(data.table)
vizDataFrame <- melt(topic_proportion_per_decade, id.vars = "year")
# plot topic proportions per decade as bar plot
library(pals)
library(ggsci)
ggplot(vizDataFrame, aes(x=year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_d3(
  )+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#Too many colors, I build my multidimensional map to distinguish between adaptation vs mitigation AND theoretical vs applied

vizDataFrame_adapt <- vizDataFrame %>% filter(variable == "V3" | variable == "V15" | variable == "V25" | variable == "V5" |
                                                variable == "V8" | variable == "V18" | variable == "V26" | variable == "V31")
ggplot(vizDataFrame_adapt, aes(x=year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_npg()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

vizDataFrame_mitigat <- vizDataFrame %>% filter(variable == "V27" | variable == "V13" | variable == "V29" | variable == "V20" |
                                                variable == "V11" | variable == "V16" | variable == "V20")
ggplot(vizDataFrame_mitigat, aes(x=year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_aaas()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

vizDataFrame_adapt_theory <- vizDataFrame %>% filter(variable == "V23" | variable == "V10" | variable == "V1" | variable == "V2" |
                                                  variable == "V9" | variable == "V12" | variable == "V22")
ggplot(vizDataFrame_adapt_theory, aes(x=year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_nejm()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

vizDataFrame_mitigat_theory <- vizDataFrame %>% filter(variable == "V19" | variable == "V7" | variable == "V4" | 
                                                         variable == "V24" | variable == "V6")
ggplot(vizDataFrame_mitigat_theory, aes(x=year, y=value, fill=variable)) + 
  geom_bar(stat = "identity") + ylab("proportion") + 
  scale_fill_jama()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


p_mod.out.corr <- topicCorr(p_stmfit2, cutoff = .01)
plot(p_mod.out.corr)

#Extract nodes and edges
p_links2 <- as.matrix(p_mod.out.corr$posadj)
p_net2 <- graph_from_adjacency_matrix(p_links2, mode = "directed")
p_net2 <- igraph::simplify(p_net2) 

p_links <- igraph::as_data_frame(p_net2, what="edges")
p_nodes <- igraph::as_data_frame(p_net2, what="vertices")

p_nodes$shape <- "dot"  
p_nodes$title <- paste0("Topic ", p_topic$TopicNumber)
p_nodes$label <- apply(p_topicNames$prob, 1, function(x) paste0(x, collapse = " \n ")) # Node label. if you want to see all the words
#Otherwise do not run 
p_nodes$size <- (p_topic$TopicProportions / max(p_topic$TopicProportions)) * 30
p_nodes$font <- "18px"
p_nodes$id <- as.numeric(1:33)
#links$color <- "gray"    # line color

visNetwork(p_nodes, p_links, width="100%",  height="800px", main="Research Topics") %>% 
  visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical")) %>%
  visNodes(scaling = list(max = 10)) %>%
  visIgraphLayout(layout ="layout_with_sugiyama", smooth = T) %>%
  visInteraction(navigationButtons = T)

visNetwork(p_nodes, p_links, width="100%",  height="800px", main="Research Topics") %>% 
  visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical")) %>%
  visNodes(scaling = list(max = 10)) %>%
  visEdges(arrows = 'from') %>%
  visIgraphLayout(physics = TRUE, smooth = T) %>%
  visHierarchicalLayout(direction = "LR", levelSeparation = 20) %>%
  visInteraction(navigationButtons = T)

plot(p_stmfit2, 
     type = "perspectives", 
     labeltype="prob",
     n = 30, 
     topics = c(17, 24), 
     text.cex = 0.8)

#Observe probabilities that each word is generated from each topic
p_td_beta <- tidy(p_stmfit2)
#Observe probabilities that each document is generated from each topic
p_td_gamma <- tidy(p_stmfit2, matrix = "gamma",
                 document_names = projects$title)

# Here I am extracting the top 10 terms per topic
p_top_terms <- p_td_beta %>%
  arrange(beta) %>%
  group_by(topic) %>%
  top_n(20, beta) %>%
  arrange(-beta) %>%
  select(topic, term) %>%
  summarise(terms = list(term)) %>%
  mutate(terms = map(terms, paste, collapse = ", ")) 

p_gamma_terms <- p_td_gamma %>%
  group_by(topic) %>%
  summarise(gamma = mean(gamma)) %>%
  arrange(desc(gamma)) %>%
  left_join(p_top_terms, by = "topic") %>%
  mutate(topic = paste0("Topic ", topic),
         topic = reorder(topic, gamma))

p_gamma_terms %>%
  top_n(20, gamma) %>%
  ggplot(aes(topic, gamma, label = terms, fill = topic)) +
  geom_col(show.legend = FALSE) +
  geom_text(hjust = 0, nudge_y = 0.0005, size = 3,
            family = "IBMPlexSans") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 0.09),
                     labels = percent_format()) +
  theme_tufte(base_family = "IBMPlexSans", ticks = FALSE) +
  theme(plot.title = element_text(size = 16,
                                  family="IBMPlexSans-Bold"),
        plot.subtitle = element_text(size = 13)) +
  labs(x = NULL, y = expression(gamma),
       title = "Top 20 topics in projects by prevalence",
       subtitle = "With the top words that contribute to each topic")

p_gamma_terms %>%
  select(topic, gamma, terms) %>%
  kable(digits = 3, col.names = c("Topic", "Expected topic proportion", "Top 7 terms"))

##Top documents per topic
p_td_gamma_top <- p_td_gamma %>%
  arrange(gamma) %>%
  group_by(topic) %>%
  top_n(8, gamma) %>%
  arrange(-gamma) %>%
  select(topic, document) %>%
  summarise(document = list(document)) %>% ##if you remove this line, you will have the top documents in decreasing order ordered one after the other
  mutate(document = map(document, paste, collapse = "; ")) %>%
  unnest(cols = document) 

final_docword <- p_td_gamma_top %>% left_join(p_top_terms, by = "topic")
kbl(final_docword) %>%
  kable_styling(bootstrap_options = c("condensed", "responsive"))

kbl(final_docword) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = TRUE, border_right = T) %>%
  column_spec(2, width = "30em") %>%
  save_kable(file = "table1.html", self_contained = T)
## Continuous effect: understand the time evolution
date <- str_split_fixed(p_meta$startDate, "/", 3)
date <- as.data.frame(apply(date, 2, as.numeric))

p_meta <- cbind(p_meta, date)

z <- projects$id[-p_processed$docs.removed]
p_stm1effect<-estimateEffect(formula=1:33 ~ frameworkProgramme,
                           stmobj=p_stmfit2,metadata=p_meta) #smooth the effect of year

p_topic.count = 33
p_model.stm.labels <- labelTopics(p_stmfit2, 1:p_topic.count)


par(mfrow=c(3,3))
for (i in seq_along(sample(1:p_topic.count, size = 33))) {
  plot(p_stm1effect, "frameworkProgramme", method = "pointestimate", topics = i, 
       main = paste0(p_model.stm.labels$prob[i,1:5], collapse = ", "), linecol = "slategrey", printlegend = F)
}

##Measure the distance between the two networks by computing the Cosine similarity between labels in nodes
library(hashr)
library(stringdist)
library(data.table)
a.int <- hash(nodes$label)
b.int <- hash(p_nodes$label)
stringmat <- stringsimmatrix(a.int, b.int, method = "cosine")
weights <- reshape2::melt(stringmat)

nodes1 <- read.csv("nodes_research.csv") #Prima era nodes_res.csv, ma non funzionava.
links1 <- read.csv("links_res.csv")
nodes_proj <- read.csv("nodes_proj.csv")
links_proj <- read.csv("links_proj.csv")

a.int <- hash(nodes1$label)
b.int <- hash(nodes_proj$label)
stringmat <- stringsimmatrix(a.int, b.int, method = "cosine")
weights <- reshape2::melt(stringmat)
weights_new <- read.csv("weights.csv")

library(ggsci)
ggplot(data = weights_new, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0.4, mid ="grey70", 
                       limits = c(0.13, +1)) +
  labs(x = "", y = "", fill = "Correlation \n Measure") +
  theme(axis.title.x = element_text(face="bold", colour="darkgreen", size = 12),
        axis.title.y = element_text(face="bold", colour="darkgreen", size = 12),
        legend.title = element_text(face="bold", colour="brown", size = 10)) +
  scale_fill_gsea()

library(heatmaply)
#The x axis is research. The y axis is project
heatmaply(stringmat, 
               dendrogram = "col",
               xlab = "", ylab = "", 
               main = "",
               scale = "column",
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = F,
               branches_lwd = 0.1,
               fontsize_row = 10, fontsize_col = 10,
               labCol = colnames(stringmat),
               labRow = rownames(stringmat),
               heatmap_layers = theme(axis.line=element_blank())
)
library(dendextend)
row_dend  <- stringmat %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 3) %>% 
  ladderize
# rotate_DendSer(ser_weight = dist(x))
col_dend  <- stringmat %>% 
  t %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram %>%
  set("branches_k_color", k = 4) %>% 
  ladderize

heatmaply_cor(
  stringmat,
  xlab = "Projects",
  ylab = "Research",
  Rowv = row_dend,
  Colv = col_dend,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "red", 
    high = "blue", 
    midpoint = 0.6, 
    limits = c(0.13, 1)
))

nodes1 <- nodes1[-3]

unif_nodes <- rbind(nodes1, nodes_proj)
names(weights_new)[1] <- "from"
names(weights_new)[2] <- "to"
names(weights_new)[3] <- "weight"
weights_new$value <- weights_new$weight
visNetwork(unif_nodes, weights_new, width="100%", height="800px", main="Distance") %>% 
  visOptions(highlightNearest = list(enabled = T, algorithm = "hierarchical")) %>%
  visNodes(scaling = list(max = 10)) %>%
  visIgraphLayout(layout ="layout_with_sugiyama", smooth = T) %>%
  visInteraction(navigationButtons = T)

corrplot(stringmat, method = "circle", cl.lim = c(0,1))
palette = colorRampPalette(c("green", "white", "red")) (20)


net_final <- graph_from_data_frame(weights_new, vertices = unif_nodes, directed = T)
c4 = cluster_optimal(net_final)
membership(c4)
coords = layout_with_fr(net_final)


