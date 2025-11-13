Netop<-function(igraph,outdir){
  # network property
  # 边数量 The size of the graph (number of edges)
  
  num.edges <- data.frame(length(E(igraph)) )# length(curve_multiple(igraph))
  num.edges
  # 顶点数量 Order (number of vertices) of a graph
  num.vertices <- data.frame(length(V(igraph)))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
  num.vertices
 
  # 连接性(connectance) 网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，可以反映网络的复杂程度
  connectance <- data.frame(edge_density(igraph,loops=FALSE))# 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
  connectance
  
  # 平均度(Average degree)
  average.degree <- data.frame((mean(igraph::degree(igraph))))# 或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
                               average.degree
                               
                               # 群连通度 edge connectivity / group adhesion
edge.connectivity <- data.frame(edge_connectivity(igraph))
edge.connectivity
                               # 聚集系数(Clustering coefficient)：分局域聚类系数和全局聚集系数，是反映网络中节点的紧密关系的参数，也称为传递性。整个网络的全局聚集系数C表征了整个网络的平均的“成簇性质”。
                               clustering.coefficient <- data.frame(transitivity(igraph))
                               clustering.coefficient
                               no.clusters <- data.frame(no.clusters(igraph))
                               no.clusters
                               # 度中心性(Degree centralization)
                               centralization.degree <- data.frame(centralization.degree(igraph)$centralization)
                               centralization.degree
                               # 介数中心性(Betweenness centralization)
                               centralization.betweenness <- data.frame(centralization.betweenness(igraph)$centralization)
                               centralization.betweenness
                               # 紧密中心性(Closeness centralization)
                               centralization.closeness <- data.frame(centralization.closeness(igraph)$centralization)
                               centralization.closeness
                               
                               num.pos.edges<-data.frame(sum(igraph::as_data_frame(igraph, what = "both")$edges$weight>0))# number of postive correlation
                               num.neg.edges<-data.frame(sum(igraph::as_data_frame(igraph, what = "both")$edges$weight<0))# number of negative correlation
                               
                               E(igraph)$weight <- abs(E(igraph)$weight)
                               
                               ############计算modularity
                               walktrap <- walktrap.community(igraph)
                               modularity <- data.frame(modularity(walktrap))
                               ############## relative_modularity 
                               relative_modularity <- data.frame(modularity / (vcount(igraph) - 1))
                               relative_modularity
                               # 平均路径长度(Average path length)
                               average.path.length <- average.path.length(igraph) # 同mean_distance(igraph) # mean_distance calculates the average path length in a graph
                               average.path.length
                               # 直径(Diameter)
                               diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
                               diameter
                               
                               
                               
                               #"average.path.length","diameter",
                               igraph.network.pro <- cbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,average.path.length,diameter,
                                                           modularity,relative_modularity,connectance,average.degree,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
                               igraph.network.pro
}