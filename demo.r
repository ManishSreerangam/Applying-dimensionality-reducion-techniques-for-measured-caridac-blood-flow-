library(dplyr)
library(ggplot2)
library(cluster) 
library(factoextra)
library(corrplot)
library(clValid)
# Preprocessing the data


data <- read.csv(file = 'data.csv', header = FALSE, row.names = 1)
head(data)
data = t(data)
ncol(data)
nrow(data)

data = as.data.frame(data)
data = dplyr::filter(data, pathology != 3)
length(data$pathology)
data.train = subset(data, select = -c(pathology))
label = data$pathology



# Viusalzing using PCA 

library(factoextra)
pcclust=prcomp(data.train,scale=FALSE) #principal component analysis
summary(pcclust)
pcclust$rotation[,1:2]

dignm<-as.character(data$pathology)
Classes = factor(label,labels = c("Healthy","BAV patients"))
pcclust$x[,1:2]


plot(pcclust$x[,1:2], col = data$pathology,xlab ="PC-1",ylab="PC-2")
legend("topright",title = "clusters",legend = unique(dignm), fill=unique(data$pathology))

# Viz using ggplot
g = ggplot(data.frame(pcclust$x[,1:2]), aes(pcclust$x[,1], pcclust$x[,2], color = Classes)) + 
  xlab("PCA 1") +
  ylab("PCA 2") +
  geom_point(size = 4)

g





# Visualizing using T-SNE plot
#perpl=35,iterations=2000

library(Rtsne)

label = as.factor(data$pathology)
Classes = factor(label,labels = c("Healthy","BAV patients"))
colors = rainbow(length(unique(label)))
names(colors) = unique(label)

tsne_plot <- function(perpl=5,iterations=500){
  set.seed(1) # for reproducibility
  tsne <- Rtsne(data.train, dims = 2, check_duplicates = FALSE,perplexity=perpl, verbose=TRUE, max_iter = iterations)
  ##Silhoutte coefficient
  df = data.frame(tsne$Y[,1],tsne$Y[,2])
  sil = silhouette(data$pathology, dist(df))
  print(mean(sil[,3]))
  ##correlation coefficient
  print(cor(tsne$Y[,1],tsne$Y[,2],  method = "pearson"))
  
  ## Connectivity
  #install.packages("clValid")
  print(connectivity(distance = dist(df), data$pathology, neighbSize = 4,
                    method = "euclidean"))
  g1 = ggplot(data.frame(tsne$Y), aes(tsne$Y[,1], tsne$Y[,2], color = Classes)) +
  xlab("TSNE 1") +
  ylab("TSNE 2")+
  ggtitle("Perplexity = 15,Iterations = 1000")+
  geom_point(size = 4)
  
}
g1 = tsne_plot()

g1

# Visualizing using UMAP
#n_neighbors = 56,min_dist = 0.5,metric = 'euclidean', random_state = as.integer(45))
#n_neighbors = 84,min_dist = 0.5,metric = 'euclidean', random_state = as.integer(12)
#n_neighbors = 12,min_dist = 0.1,metric = 'euclidean', random_state = as.integer(45)
#devtools::install_github("ropenscilabs/umapr")

library(umapr)
library(tidyverse)


df <- as.matrix(data.train)
label = as.factor(data$pathology)
Classes = factor(label,labels = c("Healthy","BAV patients"))
# run UMAP algorithm
embedding <- umap(df,n_neighbors = 56,min_dist = 0.5,random_state = as.integer(45))
#head(embedding)

g2 = embedding %>% 
  ggplot(aes(UMAP1, UMAP2, color = Classes)) +
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  ggtitle("n_neighbors = 12,min_dist = 0.1,metric = 'euclidean'")+
  geom_point(size = 4)

g2


#Silhoutte coefficient

df = data.frame(embedding$UMAP1,embedding$UMAP2)
sil = silhouette(data$pathology, dist(df))
mean(sil[,3])

#correlation coefficient
cor(embedding$UMAP1,embedding$UMAP2,  method = "pearson")

#Connectivity
#install.packages("clValid")
library(clValid)
connectivity(distance = dist(df), data$pathology, neighbSize = 4,
             method = "euclidean")

#Stability

mat = matrix(df)

stability(mat, Dist=NULL,del = 2,cluster = data$pathology, clusterDel= 2, method="euclidean")





#Working with hyperparameters

#n_neighbors
#min_dist
#n_components
#metric


neighbors <- c(4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100)



neighbors %>% 
  map_df(~umap(as.matrix(data.train), n_neighbors = .x, random_state = as.integer(45)) %>% 
           mutate(pathology = data$pathology, Neighbor = .x)) %>% 
  mutate(Neighbor = as.integer(Neighbor)) %>% 
  ggplot(aes(UMAP1, UMAP2, color = pathology)) + 
  geom_point() + 
  facet_wrap(~ Neighbor, scales = "free")







# Multi Dimensional Scaling

data.train <- dist(data.train) # euclidean distances between the rows
fit <- cmdscale(data.train,eig=TRUE, k=2, add = FALSE, x.ret = FALSE)

x <- fit$points[,1]
y <- fit$points[,2]

ggplot(data.frame(fit$points), aes(x, y, color = Classes)) + geom_point(size = 3)


g1


install.packages("ggpubr")
library(ggpubr)

ggarrange(g,g1,g2,g3, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)


library(ggplot2)
library(Rtsne)
tsne_plot <- function(perpl=30,iterations=2000){
  set.seed(1) # for reproducibility
  tsne <- Rtsne(data.train, dims = 2, check_duplicates = FALSE,perplexity=perpl, verbose=TRUE, max_iter = iterations)
  ggplot(data.frame(tsne$Y), aes(tsne$Y[,1], tsne$Y[,2], color = data$pathology)) + geom_point()
  
}

perplexity_values <- c(2,5,30,50,100)
sapply(perplexity_values,function(i){tsne_plot(perpl=i)})

#install.packages("devtools")
library(usethis)
library(devtools)
install_github("crj32/M3C")
library(tsne)
library(M3C)


result = tsne(df, data$pathology, perplex=5)
result = data.frame(result)

ggplot(result, aes(X1, X2, color = data$pathology)) + geom_point()


p <- ggplot(mtcars, aes(mpg, wt)) +
  geom_point(aes(colour = factor(cyl)))
p + scale_colour_manual(values = c("red", "blue", "green"))




