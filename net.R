#################network


otu<- read.csv("otutable.csv",row.names = 1)
dim(otu)
otu_rare<-otu

trt_id <- c('SoilSoil','RhizosphereHealthy','RhizosphereDisease',
            'RootHealthy','RootDisease','StemHealthy','StemDisease',
            'SeedHealthy','SeedDisease')
col_names <- sapply(trt_id, function(x) grep(x, colnames(otu_rare)))
split_otu <- lapply(col_names, function(x) {
  x_data <- otu_rare[,x]
  x_data <- x_data[-which(rowSums(x_data)==0),]
})

tax<-read.csv("futax.csv",row.names = 1)

g <- lapply(split_otu,function(x){###加权
  x<-data.frame(x)
  x<-x[which(rowMeans(x) > 0.0001), ]
  x<-cbind(x,tax[row.names(x),])
  x1<-x[1:c(ncol(x)-7)]
  b<-rcorr(t(x1),type="spearman")
  rr<-b$r
  pp<-b$P
  rr[abs(rr)<0.4]<-0
  pp<-p.adjust(pp,method="BH")
  pp[pp>=0.05&pp<1]<-0  
  pp[pp<0.05]<-1
  z<-rr*pp
  diag(z)<-0
  g<-graph.adjacency(z,weighted=TRUE,mode="undirected")
  # 删除自相关
  g <- simplify(g)
  p<-x[,ncol(x)-5]
  c<-x[,c(ncol(x)-4)]
  o<-x[,c(ncol(x)-3)]
  f<-x[,c(ncol(x)-2)]
  ge<-x[,c(ncol(x)-1)]
  V(g)$phylum<-p
  V(g)$class<-c
  V(g)$order<-o
  V(g)$family<-f
  V(g)$genus<-ge
  # 删除孤立节点
  g <- delete.vertices(g, which(degree(g)==0) )
  V(g)$degree<-igraph::degree(g)
  V(g)$Cen_degree<-centralization.degree(g)$res
  V(g)$Ben_degree<-centralization.betweenness(g)$res
  V(g)$Clo_degree<-centralization.closeness(g)$res
  return(g)
})

save(g,file = 'network0.0001r0.4p0.05.rda')

library(igraph)
library(dplyr)
library(Hmisc)
library(WGCNA)
library(preprocessCore)

#library(rsparcc)

load('network0.0001r0.4p0.05.rda')

cols<-c("#f49128","#194a55","#187c65","#f26115","#c29f62","#83ba9e","#c62d17",
        "#023f75","#ea894e","#266b69","#eb4601","#f6c619","#fa6e01","#2f2f2f",
        "#972b1d","#e6a84b","#4c211b","#ff717f")
col_g <- "#C1C1C1"

pdf(paste0("ba-pre.pdf"), encoding="MacRoman", width=15, height=9)
par(mfrow=c(2,5),mar=c(0,0,1,0),font.main=4)
for(i in 1:9){
  g1 <- g[[i]]
  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)
  set.seed(007)
  V(g1)$modularity <- membership(cluster_fast_greedy(g1))
  
  V(g1)$label <- V(g1)$name
  V(g1)$label <- NA
  modu_sort <- V(g1)$modularity %>% table() %>% sort(decreasing = T)
  
  top_num <- 8
  modu_name <- names(modu_sort[1:8])
  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  V(g1)$color <- V(g1)$modularity
  V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
  V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)],modu_name)]
  V(g1)$frame.color <- V(g1)$color
  
  E(g1)$color <- col_g
  for ( i in modu_name){
    col_edge <- cols[which(modu_name==i)]
    otu_same_modu <-V(g1)$name[which(V(g1)$modularity==i)]
    E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
  }
  
  sub_net_layout <- layout_with_fr(g1, niter=999,grid = 'nogrid')
  plot(g1,layout=sub_net_layout, edge.color = E(g1)$color,vertex.size=2)
  
  title(main = paste0('Nodes=',length(V(g1)$name),', ','\nEdges=',nrow(data.frame(as_edgelist(g1)))))
  legend("right", legend = modu_name, col = modu_cols, bty = "n", fill = modu_cols, border = modu_cols)  # Update legend parameters
}
dev.off()


##################Fig1A-B


dis<-read.csv("Dis.csv")
rhd <- dis %>% 
  filter(Compartment == "Rhizosphere")
rhd <-rhd %>% 
  filter(State == "Disease")
row.names(rhd)<-rhd[,1]
coord<-rhd[6:7]
env<-rhd[9:22]

require(ape)
tree1<-read.tree("rhd_tree.nwk") 
comstate2<-list()
for (i in 1:length(rownames(positive_pairs))){
  minicom0<-minicom[rownames(minicom)==positive_pairs[i,1]|rownames(minicom)==positive_pairs[i,2],]
  comstate<-matrix(NA, ncol=4, nrow=length(colnames(minicom)))
  
  for (j in 1:length(colnames(minicom)))
  {
    comstate[j,]<-c(paste(as.character(positive_pairs[i,1]), as.character(positive_pairs[i,2]), sep="_"), paste(as.numeric(positive_pairs[i,3]), sep=""), colnames(minicom)[j], paste(as.vector(minicom0[,j])[1], as.vector(minicom0[,j])[2], sep=""))
  }
  comstate2[[i]]<-comstate
}

result<-do.call(rbind.data.frame, comstate2)
colnames(result)<-c("pair","weight","com","state")
isotopic<-result[which(result$state=="00"|result$state=="11"), ] 



pca_env <-prcomp(env, scale.=TRUE) 
env$PC1_env <- pca_env$x[,1]
env$PC2_env <- pca_env$x[,2] 

env$com<-rownames(env)
coord$com<-rownames(coord)
megatable_pos<-merge(isotopic, coord, by="com") 
megatable_pos<-merge(megatable_pos, env, by="com") 

megatable_pos1<-as.data.frame(megatable_pos)



megatable_pos$pair<-factor(megatable_pos$pair)



pvalue<-NULL
for (i in 1:length(levels(megatable_pos$pair)))
{
  dataset<-megatable_pos[which(megatable_pos$pair %in% levels(megatable_pos$pair)[i]),]
  ifelse(length(levels(factor(dataset$state)))==1, pvalue[i]<-1,
         {
           Y<-with(dataset, cbind(UTM_X, UTM_Y)) #Use spatial coordinates variables
           
           mymanova<-manova(Y~dataset$state) # For each aggregated pair, is the geographic distance significantly different between communities with isotopic states?  
           sum_mymanova<-unlist(summary(mymanova))
           pvalue[i]<-as.numeric(sum_mymanova["stats11"]) #Extract the disp_pvalue
           
         })
}

# Final table for aggregated pairs
final_pos<-cbind.data.frame(as.character(positive_pairs$sp1), as.character(positive_pairs$sp2), as.numeric(round(positive_pairs$weight,4)), round(pvalue, 4)) # Add the disp_pvalue to "final_pos"
colnames(final_pos)<-c("sp1", "sp2", "weight", "disp_pvalue")
final_pos$process1<-ifelse(final_pos$disp_pvalue<0.05, "dispersal_limitation", "no_dispersal_limitation")

# 3. Testing the existence of environmental filtering vs. positive interaction 
##############################################################################

# Perform MANOVAs across all the pairs of species involved in aggregated associations

pvalue<-NULL
for (i in 1:length(levels(megatable_pos$pair)))
{
  dataset<-megatable_pos[which(megatable_pos$pair %in% levels(megatable_pos$pair)[i]),]
  ifelse(length(levels(factor(dataset$state)))==1, pvalue[i]<-1,
         {
           Y<-with(dataset, cbind(PC1_env, PC2_env)) #Use environmental variables
           
           mymanova<-manova(Y~dataset$state) # For each aggregated pair, is the environmental data significantly different between communities with isotopic states? 
           sum_mymanova<-unlist(summary(mymanova))
           pvalue[i]<-as.numeric(sum_mymanova["stats11"]) #Extract the filtering_pvalue
           
         })
}

# Final table for aggregated pairs
final_pos_env<-cbind.data.frame(as.character(positive_pairs$sp1), as.character(positive_pairs$sp2), round(pvalue, 4)) # Add the filtering_pvalue to "final_pos_env"
colnames(final_pos_env)<-c("sp1", "sp2", "filtering_pvalue")
final_pos_env$process2<-ifelse(final_pos_env$filtering_pvalue<0.05, "filtering", "positive_int")
final_pos$filtering_pvalue<-final_pos_env$filtering_pvalue 

#Add the filtering_pvalue to "final_pos"
final_pos$process2<-final_pos_env$process2 
final_pos$final_process<-with(final_pos, ifelse(process1=="dispersal_limitation", ifelse(process2=="filtering", "Disp_OR_filtering", "dispersal_limitation"), process2)) # If disp_pvalue and filtering_pvalue are <0.06, both dispersal limitation and environmental filtering could explain the aggregated pair


# 4. Adding phylogenetic distances to the final table
#######################################################

require(ape)
cophenetic1<-cophenetic(tree1) #Calculate pairwise phylogenetic distances

phydist1<-NULL
for (i in 1:nrow(final_pos))
{
  phydist1[i]<-cophenetic1[as.character(final_pos$sp1[i]), as.character(final_pos$sp2[i])]
}
final_pos<-cbind(final_pos, phydist1) # Add phydist1 (in Myr) to "final_pos"

# Repeat with replicated trees if needed

write.csv(final_pos,"????Assembly_AggregatedPairs.csv")




comstate2<-list()

for (i in 1:length(rownames(negative_pairs)))
{
  minicom0<-minicom[rownames(minicom)==negative_pairs[i,1]|rownames(minicom)==negative_pairs[i,2],]
  comstate<-matrix(NA, ncol=4, nrow=length(colnames(minicom)))
  
  for (j in 1:length(colnames(minicom)))
  {
    comstate[j,]<-c(paste(as.character(negative_pairs[i,1]), as.character(negative_pairs[i,2]), sep="_"), paste(as.numeric(negative_pairs[i,3]), sep=""), colnames(minicom)[j], paste(as.vector(minicom0[,j])[1], as.vector(minicom0[,j])[2], sep=""))
  }
  comstate2[[i]]<-comstate
}

result<-do.call(rbind.data.frame, comstate2)
colnames(result)<-c("pair","weight","com","state")
allotopic<-result[which(result$state=="01"|result$state=="10"), ] # Relevant for segregated pairs


# Megatable for segregated pairs
megatable_neg<-merge(allotopic, coord, by="com") 
megatable_neg<-merge(megatable_neg, env, by="com") # Merge "allotopic" with "coord" and "env"


# 2. Testing the existence of dispersal limitation 
###################################################
megatable_neg$pair<-factor(megatable_neg$pair)
# Perform MANOVAs across all the pairs of species involved in segregated associations

pvalue<-NULL
for (i in 1:length(levels(megatable_neg$pair)))
{
  dataset<-megatable_neg[which(megatable_neg$pair %in% levels(megatable_neg$pair)[i]),]
  ifelse(length(levels(factor(dataset$state)))==1, pvalue[i]<-1,
         {
           Y<-with(dataset, cbind(UTM_X, UTM_Y)) #Use spatial coordinates variables
           
           mymanova<-manova(Y~dataset$state) # For each segregated pair, is the geographic distance significantly different between communities with allotopic states?  
           sum_mymanova<-unlist(summary(mymanova))
           pvalue[i]<-as.numeric(sum_mymanova["stats11"]) #Extract the disp_pvalue
           
         })
}

# Final table for segregated pairs
final_neg<-cbind.data.frame(as.character(negative_pairs$sp1), as.character(negative_pairs$sp2), as.numeric(round(negative_pairs$weight,4)), round(pvalue, 4)) # Add the disp_pvalue to "final_neg"
colnames(final_neg)<-c("sp1", "sp2", "weight", "disp_pvalue")
final_neg$process1<-ifelse(final_neg$disp_pvalue<0.05, "dispersal_limitation", "no_dispersal_limitation")


# 2. Testing the existence of environmental filtering vs. negative interaction 
##############################################################################

# Perform MANOVAs across all the pairs of species involved in segregated associations

pvalue<-NULL
for (i in 1:length(levels(megatable_neg$pair)))
{
  dataset<-megatable_neg[which(megatable_neg$pair %in% levels(megatable_neg$pair)[i]),]
  ifelse(length(levels(factor(dataset$state)))==1, pvalue[i]<-1,
         {
           Y<-with(dataset, cbind(PC1_env, PC2_env)) #Use environmental variables
           
           mymanova<-manova(Y~dataset$state) # For each segregated pair, is the environmental data significantly different between communities with allotopic states? 
           sum_mymanova<-unlist(summary(mymanova))
           pvalue[i]<-as.numeric(sum_mymanova["stats11"]) #Extract the filtering_pvalue
           
         })
}

# Final table for segregated pairs
final_neg_env<-cbind.data.frame(as.character(negative_pairs$sp1), as.character(negative_pairs$sp2), round(pvalue, 4)) # Add the filtering_pvalue to "final_neg_env"
colnames(final_neg_env)<-c("sp1", "sp2", "filtering_pvalue")
final_neg_env$process2<-ifelse(final_neg_env$filtering_pvalue<0.05, "filtering", "negative_int")
final_neg$filtering_pvalue<-final_neg_env$filtering_pvalue 


#Add the filtering_pvalue to "final_neg"
final_neg$process2<-final_neg_env$process2 
final_neg$final_process<-with(final_neg, ifelse(process1=="dispersal_limitation", ifelse(process2=="filtering", "Disp_OR_filtering", "dispersal_limitation"), process2)) # If disp_pvalue and filtering_pvalue are <0.06, both dispersal limitation and environmental filtering could explain the segregated pair


# 4. Adding phylogenetic distances to the final table
######################################################

phydist1<-NULL
for (i in 1:nrow(final_neg))
{
  phydist1[i]<-cophenetic1[as.character(final_neg$sp1[i]), as.character(final_neg$sp2[i])]
}
final_neg<-cbind(final_neg, phydist1) # Add phydist1 (here, in Myr) to "final_neg"

write.table(final_neg, file="Assembly_SegregatedPairs.csv")


##############################################################
# Testing statistical differences in pairwise phylogenetic   # 
# relatedness across the four main assembly processes        #
##############################################################


# Combining final_pos and final_neg to test phydist differences across the four main assembly processes
final_pos2<-final_pos
final_neg2<-final_neg
final_pos2$final_process<-ifelse(final_pos2$final_process=="filtering", "positive_filtering", final_pos2$final_process)
final_neg2$final_process<-ifelse(final_neg2$final_process=="filtering", "negative_filtering", final_neg2$final_process)
final_tot<-rbind(final_pos2, final_neg2)
final_tot<-final_tot[-which(final_tot$final_process=="Disp_OR_filtering"),]
final_tot<-final_tot[-which(final_tot$final_process=="dispersal_limitation"),]


# Run permutation anova and pairwise permutation t-tests
#install.packages("RVAideMemoire")
require(RVAideMemoire)
perm.anova(final_tot$phydist1~final_tot$final_process)
pairwise.perm.t.test(final_tot$phydist1, final_tot$final_process)

head(final_tot)
write.csv(final_tot,"final_tot.csv")




