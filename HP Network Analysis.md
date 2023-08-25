---
title: "Harry Potter Network Analysis"
author: "Federico Moroni"
date: "2023-07-06"
output:
  md_document:
    variant: markdown_github
---
The **Harry Potter Network Analysis** aims to delve into the intricate web of characters and connections within the Harry Potter universe.Beyond the surface level of this beloved fantasy series lies a complex web of relationships, interactions, and influences that extend far beyond the pages of the books or the frames of the films.

This study aims to embark on a comprehensive network analysis of the Harry Potter universe, exploring the complex network of characters and their interactions. By employing sophisticated network analysis techniques, we seek to unveil the underlying structure and dynamics of this fantastical realm, shedding new light on the intricate relationships between its various components.

At the heart of this analysis lies the concept of network theory, which provides a powerful framework for understanding complex systems. By representing the Harry Potter universe as a network, where individual elements are nodes and their relationships are represented as edges, we can unravel the intricate network of connections that give life to this fictional world.

---

<H1 align=center>Structure</H1>

0. *Network Representation*: transform the data into a suitable network representation. Represent the entities as nodes and the relationships as edges or links between the nodes. 

1. *Network Visualization*: our first objective is to visually represent the intricate web of connections within the Harry Potter Network. By creating a visual representation, we can gain a comprehensive overview of the character connections, identify key figures, and discover clusters or groups of characters with significant relationships.

2. *Network Analysis*: through the calculation of centrality measures, such as degree centrality, betweenness centrality, closeness centrality, and pagerank centrality, we aim to identify the most influential and prominent characters in the Harry Potter Network. This analysis will reveal characters who play crucial roles in connecting different parts of the network and those who have the highest overall influence within the narrative.

3. *Structural Analysis*: we will explore the structural properties of the Harry Potter Network, utilizing statistics such as diameter, mean distance and transitivity. By examining these measures, we can describe the overall cohesion and complexity of the network, uncover recurring patterns, and reveal narrative themes and character dynamics.

4. *Network modeling*: to understand the underlying processes that generate the connections within the Harry Potter Network, we will employ an Exponential Random Graph Model (ERGM). By fitting the model to our dataset, we can estimate the probabilities of specific types of connections and gain insights into the mechanisms that shape the network's structure.

---

```{r include=FALSE}
library(igraph)
library(ape)
library(dplyr)
library(visNetwork)
library(ergm)
library(dplyr)
```

```{r}
nodes = read.csv('characters.csv')
table(is.na(nodes))
```

```
## 
## FALSE 
##   195
```

```{r}
edges = read.csv('relations.csv')
table(is.na(edges))
```

```
## 
## FALSE 
##   1539
```

It is crucial to examine the presence of missing values within our dataset before commencing the analysis. 

Fortunately, in our specific case, we do not encounter any missing values, ensuring the completeness and reliability of our data for analysis.

<H2 align=center>0-NETWORK REPRESENTATION</H2> 

```{r}
#We build an undirect graph object
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
```

```{r}
paste("The network has",ecount(g),"edges")
```
```
## [1] "The network has 513 edges"
```
```{r}
paste("The network has",vcount(g),"nodes")
```
```
## [1] "The network has 65 nodes"
```
```{r}
#Check if the graph is connected
is.connected(g)
```
```
## [1] TRUE
```
A connected graph is characterized by the presence of a path between every pair of vertices. In simpler terms, in a connected graph, it is possible to reach any vertex from any other vertex by traversing a sequence of edges.

```{r}
#Check if the edges of the graph are weighted
is.weighted(g)
```
```
## [1] FALSE
```
```{r}
#Check the simplicity of the graph
is.simple(g)
```
```
## [1] FALSE
```
We can see that the graph is not simple, hence we can have self-loops or multiple edges between the same pair of vertices.

<H2 align=center>1-NETWORK VISUALIZATION</H2>

```{r}
#Plot the graph
igraph_options(vertex.size=4,edge.arrow.size = 0.5)
plot(g, layout = layout_with_lgl, vertex.label = "")
```
![unnamed-chunk-10-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/4e921de4-ac2c-4553-8b54-4fd87ca6a5b8)


<H2 align=center>2-NETWORK ANALYSIS</H2>

Since there are various methods to determine centrality which can yield to distinct results, we will draw comparison between the diffent methodologies.

### a. Degree distribution

The degree of a node in a graph corresponds to the number of edges connecting to that specific vertex. Essentially, it indicates the count of direct connections a particular vertex has with neighboring vertices. Nodes with a high degree tend to be more influential in the network and may play a key role in information
diffusion or interaction among nodes.

```{r}
D <- igraph::degree(g)
sort(D, decreasing = TRUE)[1:10]
```
```
## Hermione Granger      Ron Weasley     Harry Potter Albus Dumbledore 
##               59               57               56               55 
##   Lord Voldemort   George Weasley   Arthur Weasley     Sirius Black 
##               46               35               34               33 
##    Ginny Weasley     Bill Weasley 
##               30               28
```
By examining the degrees of our nodes we observe that **Hermione Granger** possesses the highest degree, hence she is the most connected node within our graph.

```{r}
plot(sort(D), main = "Degree Distribution", ylab = "Degree", col="blue")
```

![unnamed-chunk-12-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/e7db0208-d30d-4e95-83ee-c07aa18a766b)

Particularly we can identify a group of four characters (**Hermione Granger**, **Ron Weasley**, **Harry Potter** and **Albus Dumbledore**) that are the individuals who hold most information or can quickly connect with the wider network.

```{r}
#Plot the histogram of the distribution of degree
hist(D, xlab="Degree", main="Histogram of degree distribution", col = "blue")
#Add the mean
abline(v = mean(D), col = "red", lwd=2)
#Add the median
abline(v = median(D), col = "green", lwd=2)
```
![unnamed-chunk-13-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/2e74a8ee-d812-4937-9948-b648f6a1885d)

From the histogram we can see that degree distribution is very concentrated in the first interval (0-10), in particular on the left of the mean degree.
Comparing the mean (read line) and the median (green line), it is possible to confirm that the distribution is not symmetric.
It is evident that we have a lot of nodes with small degree and a few nodes with very large degree. This suggest that the degree distribution is heterogeneous.

### b. Betweeness

Betweenness centrality measures the extent to which a node serves as a bridge or intermediary between other nodes. It quantifies the number of shortest paths that pass through a node, highlighting its role in facilitating communication and information flow within the network. Nodes with high betweenness centrality are crucial connectors, controlling the flow of information and exerting influence over the network.

```{r}
betw <- betweenness(g)
sort(betw, decreasing = TRUE)[1:10]
```
```
##        Harry Potter      Lord Voldemort         Ron Weasley    Albus Dumbledore 
##           481.93569           386.41922           238.86777           225.89852 
##    Hermione Granger  Neville Longbottom       Luna Lovegood        Sirius Black 
##           206.19561            83.00306            63.00000            61.45676 
##       Rubeus Hagrid Bellatrix Lestrange 
##            44.59529            31.31766
```
Betweenness indicates that **Harry Potter** is the one that acts as a critical intermediary or bridge.

```{r}
plot(sort(betw), main = "Betweenness", ylab = "Betweenness", col="blue")
```
![unnamed-chunk-15-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/2f331075-0410-4de9-ac59-90c5095b30bc)

```{r}
hist(betw, col="blue", main="Histogram of Betweenness", xlab="Betweenness")
abline(v = mean(betw), col = "red", lwd=2)
abline(v = median(betw), col = "green", lwd=2)
```
![unnamed-chunk-16-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/5154e538-6d71-4ccd-b1c7-8eebe9d1493a)

From the plots we can clearly see that there are very few nodes with high centrality, hence Harry Potter is by far the most influential character.
Instead, apart from a few exceptions, the other characters have very little influence.

### c. Closeness

Closeness centrality quantifies how closely connected a node is to other nodes in a network. It is calculated as the inverse of the average distance from a node to all other nodes in the network. Nodes with high closeness centrality are generally more accessible and have shorter average path lengths to other nodes. They are well-positioned to efficiently spread information or influence throughout the network. Closeness centrality helps identify nodes that can quickly disseminate or receive information and play a significant role in information diffusion and network efficiency.

```{r}
close <- closeness(g, mode="total")
sort(close, decreasing = TRUE)[1:10]
```
```
##     Harry Potter      Ron Weasley Hermione Granger   Lord Voldemort 
##      0.012500000      0.010989011      0.010869565      0.010752688 
## Albus Dumbledore     Sirius Black    Ginny Weasley     Fred Weasley 
##      0.010204082      0.009259259      0.008928571      0.008849558 
##   George Weasley    Molly Weasley 
##      0.008849558      0.008849558
```
We can see that **Harry Potter** has the highest value, hence, looking at closeness, he is the central node.

```{r}
plot(sort(close), main = "Closeness", ylab = "Closeness", col="blue")
```
![unnamed-chunk-18-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/0f718483-424f-4719-b016-2a6094ff27ee)

```{r}
hist(close, col = "blue", main="Histogram of Closeness", xlab="Closeness")
abline(v = mean(close), col = "red", lwd=2)
abline(v = median(close), col = "green", lwd=2)

```
![unnamed-chunk-19-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/b0f3d12a-4a5d-4277-839b-7a93bdffb9c4)

From the plots we can clearly see that another time Harry Potter is by far the most connected character.

### d. Pagerank

Pagerank is used to measure the importance of nodes based on their connections. It assigns scores to nodes based on the number and quality of incoming links, with highly connected and influential nodes receiving higher scores. Pagerank helps identify important nodes in large networks and is useful for understanding influence and information flow within a network.

```{r}
prank <- page_rank(g)
sort(prank$vector, decreasing = TRUE)[1:10]
```
```
##     Harry Potter Hermione Granger      Ron Weasley Albus Dumbledore 
##       0.05414539       0.05093100       0.04992041       0.04902116 
##   Lord Voldemort     Sirius Black   George Weasley   Arthur Weasley 
##       0.04375762       0.02838925       0.02676012       0.02615194 
##    Ginny Weasley     Fred Weasley 
##       0.02378409       0.02204426
```
From Pagerank we can see that **Harrry Potter** is the most important character.

```{r}
plot(sort(prank$vector), main = "Pagerank", ylab = "Pagerank", col="blue")
```
![unnamed-chunk-21-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/71233739-7eee-460d-8349-cd15b9478e00)

```{r}
hist(prank$vector, col="blue", main="Histogram of Pagernk", xlab="Pagerank")
abline(v = mean(prank$vector), col = "red", lwd=2)
abline(v = median(prank$vector), col = "green", lwd=2)
```
![unnamed-chunk-22-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/400a61ad-f100-4729-9eea-0a49d59e4d44)

From the graphs we can see that there is a small group of important characters.

<H3 align=center>2.1-Power Law Analysis</H3>

```{r}
degree_dist <- function (g) {
 fd <- table(igraph::degree(g))
 d <- as.numeric(names(fd)) + 1 # degree + 1
 list(d = d, fd = fd)
}
```

```{r}
fd <- degree_distribution(g)
d <- igraph::degree(g)
dd <- degree_dist(g)

(m0 <- lm(log(fd) ~ log(d) , data = dd))
```
```
## 
## Call:
## lm(formula = log(fd) ~ log(d), data = dd)
## 
## Coefficients:
## (Intercept)       log(d)  
##      1.8398      -0.4752
```
```{r}
(m1 <- glm(fd ~ log(d), family = poisson, data = dd))
```
```
## 
## Call:  glm(formula = fd ~ log(d), family = poisson, data = dd)
## 
## Coefficients:
## (Intercept)       log(d)  
##      2.0489      -0.4996  
## 
## Degrees of Freedom: 31 Total (i.e. Null);  30 Residual
## Null Deviance:       31.85 
## Residual Deviance: 19.46     AIC: 101.4
```
```{r}
with(dd, plot(log(d), log(fd), main="Linear regr. vs Poisson log-linear regr."))
abline(m0$coefficients[1], m0$coefficients[2], col = 'red')
abline(m1$coefficients[1], m1$coefficients[2], col = 'blue')
```
![unnamed-chunk-26-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/975dff3d-35ef-4bea-910a-79ef8192be4d)

We plot two power law exponent estimates, one from linear regression and another from Poisson log-linear regression.

The two lines are clearly parallel; the blue line (Poisson model) seems to fit better, but there are no great differences between the values of intercept and slope.

<H2 align=center>3-STRUCTURAL ANALYSIS</H2>

### a. Diameter

Diameter represents the maximum shortest path length between any two nodes in a network. It measures the longest distance or the maximum number of edges one must traverse to travel between the farthest pair of nodes in the network. The diameter provides insights into the overall size and connectivity of a network, indicating the maximum potential distance for communication or information transmission within the network.

```{r}
diam <- diameter(g)
diam
```
```
## [1] 4
```
The maximum shortest path length between any two nodes in the network is relatively short. This indicates a high level of connectivity and efficient communication within the network. It implies that information or influence can spread quickly across the network, allowing for efficient interactions and rapid dissemination of information.

### b. Mean Distance

Mean distance refers to the average shortest path length between all pairs of nodes in a network. It provides a measure of the overall average distance or number of edges required to travel between any two nodes in the network. The mean distance offers insights into the average proximity and connectivity within the network, indicating how closely or distantly connected the nodes are on average. It helps assess the efficiency of information transmission, communication, and overall network accessibility.

```{r}
mean_distance(g)
```
```
## [1] 2.028365
```
A mean distance of 2 indicates that, on average, nodes in the network are very closely connected. It suggests a highly efficient and well-connected network structure, where most nodes can be reached within just a few steps. Such a low mean distance implies quick and direct communication, efficient information flow, and easy accessibility between nodes

### c. Clustering

Clustering refers to the tendency of nodes in a network to form densely connected groups or clusters. It measures the degree to which nodes within a cluster are more strongly connected to each other than to nodes in other clusters. Clustering helps identify cohesive subgroups within a network, where nodes have a higher density of connections within their own group.

```{r}
transitivity(g)
```
```
## [1] 0.4133759
```
```{r}
hist(transitivity(g, "local"), col="blue", main="Histogram of Transitivity", xlab="Transitivity")
```
![unnamed-chunk-30-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/bd1bf1a5-b83f-4d12-b527-71be0d29ae85)

Transitivity reflects the presence of clustering or triadic closure in a network. It suggests that if two nodes in the network share a common connection, they are likely to have additional connections with each other. This phenomenon can be seen as a measure of the overall connectedness and cooperative relationships within the network.
Generally, transitivity values range from 0 to 1, with higher values indicating a greater degree of transitivity or clustering in the network.

In this case we have a value of 0.41 that suggests a moderate level of clustering or transitivity in the network, indicating that there are some triangles or closed triplets present, but the network is not fully transitive.

<H3 align=center>3.1-GRAPH PARTITIONING</H3>

### a. Community detection based on greedy optimization of modularity

Method used to identify cohesive and distinct communities within a network.
Modularity is a measure that quantifies the strength of the division of a network into communities, by comparing the density of connections within communities to the expected density in a random network.
The goal of this method is to find a partitioning of the network that maximizes the modularity, indicating the presence of meaningful and well-defined communities.

```{r}
#Having a non-simple graph, firstly we have to simplify it
g <- simplify(g)
kc <- fastgreedy.community(g)
sizes(kc)
```
```
## Community sizes
##  1  2  3  4 
## 15 27  8 15
```
We have 4 groups with similar sizes.

```{r}
fc <- cluster_fast_greedy(g)
head(membership(fc))
```
```
## Regulus Arcturus Black           Sirius Black         Lavender Brown 
##                      1                      2                      4 
##              Cho Chang     Vincent Crabbe Sr.         Vincent Crabbe 
##                      4                      1                      1
```
Looking at this output, it is possible to see the community of each hero.

```{r}
plot(kc, g, vertex.label = "",asp=0.5,vertex.size=3)
```
![unnamed-chunk-33-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/e0195029-28c7-4d67-9606-efc28995e66f)

```{r}
plot_dendrogram(kc, mode= "phylo")
```
![unnamed-chunk-34-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/530ce7f4-5b50-4b1d-acb6-5532a68db96c)

It is based on the modularity measure and a hierarchical approach.

### b. Louvain Comunity Detection

The Louvain Community Detection algorithm is based on optimizing modularity, which measures the strength of division of a network into communities. The algorithm iteratively optimizes modularity by merging and rearranging nodes between communities to maximize the modularity score.
It has become widely adopted due to its computational efficiency and ability to identify communities with high modularity. However, it is important to note that, like other community detection algorithms, is a heuristic method and does not guarantee finding the globally optimal community structure.

```{r}
nodes <- nodes %>%
  rename("id" = "id",
         "label" = "name")

edges <- edges %>% 
        rename("from" = "source",
               "to" = "target")
```


```{r}
cluster <- cluster_louvain(g)

cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- nodes$label

# Create group column
nodes <- left_join(nodes, cluster_df, by = "label")
colnames(nodes)[4] <- "group"
```

```{r}
visNetwork(nodes, edges, width = "100%") %>%
  visIgraphLayout() %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#0085AF",
      border = "#013848",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = TRUE, size = 10)
  ) %>%
  visEdges(
    shadow = FALSE,
    color = list(color = "#0085AF", highlight = "#C62F4B")
  ) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
             selectedBy = "group") %>% 
  visLayout(randomSeed = 11)
```

### c. Main Characters Interactions

```{r}
# who is the last
ia <- order(betw, decreasing = TRUE)[1]
V(g)$name[ia]
```
```
## [1] "Harry Potter"
```
```{r}
# give us the position
tail(sort(betw))
```
```
## Neville Longbottom   Hermione Granger   Albus Dumbledore        Ron Weasley 
##           83.00306          206.19561          225.89852          238.86777 
##     Lord Voldemort       Harry Potter 
##          386.41922          481.93569
```
```{r}
# ego network
sg <- subgraph.edges(g, E(g)[.inc(ia)][1:10])
plot(sg)
```
![unnamed-chunk-40-1](https://github.com/federicomoroni/Harry-Potter-Network-Analysis/assets/102738225/e6d945e8-9cb8-4af5-b8cd-e8b36d2c2a40)

The ego network provides a localized view of the relationships and interactions that surround the ego (central node) within the larger network.

<H2 align=center>4-NETWORK MODELING</H2>

The Exponential Random Graph Model is a statistical model used to analyze and explain the patterns of connections in a network.
It is based on the principles of exponential families in statistics and use a set of parameters to model the probability of observing a specific network configuration

```{r message=FALSE}
am <- get.adjacency(g, sparse = FALSE)
g1 <- as.network(am, directed = FALSE)
ergm(g1~edges) %>% summary()
```
```
## Call:
## ergm(formula = g1 ~ edges)
## 
## Maximum Likelihood Results:
## 
##       Estimate Std. Error MCMC % z value Pr(>|z|)    
## edges -1.66828    0.06001      0   -27.8   <1e-04 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##      Null Deviance: 2883  on 2080  degrees of freedom
##  Residual Deviance: 1820  on 2079  degrees of freedom
##  
## AIC: 1822  BIC: 1827  (Smaller is better. MC Std. Err. = 0)
```
- The negative estimate for the edges coefficient suggests that the presence of edges in the network has a strong negative effect on the likelihood of the observed network structure. This means that as the number of edges in the network increases, the probability of observing the specific pattern or configuration of edges in the analyzed network decreases.

- The MCMC percentage indicates the proportion of MCMC samples used for model estimation. In this case, it is 0, indicating that the estimation was done using maximum likelihood estimation (MLE).

- The p-value is highly significant, hence the negative relationship between the presence of edges and the network structure is significant.

- the reduction in deviance indicates that the fitted model provides a better fit to the data compared to the null model.

Overall, the output suggests that the ‘edges’ term in the ERGM model has a significant negative effect on the network. The estimate indicates a strong negative relationship between the presence of edges and the network structure. The small p-value suggests that this effect is unlikely to occur by chance.

<H1 align=center>Conclusion</H1>

In conclusion, our network analysis of the Harry Potter universe has provided valuable insights into the intricate web of connections that underlies this captivating fictional world. Through network visualization, we gained an overview of character connections and identified key figures and clusters within the network.

By calculating centrality measures, such as degree centrality and closeness centrality, we identified influential and prominent characters, uncovering their significant roles within the Harry Potter Network.

Furthermore, the application of an Exponential Random Graph Model (ERGM) helped us understand the underlying generating processes that shape the network. By fitting the model to our dataset, we gained insights into the probabilities of specific types of connections, shedding light on the network's formation and dynamics.

This network analysis has deepened our understanding of the Harry Potter universe, revealing the complex relationships, character dynamics, and narrative arcs that make this fictional world so captivating. The insights gained from this analysis can contribute to the appreciation and exploration of the rich tapestry of characters and their interactions within the Harry Potter series.

