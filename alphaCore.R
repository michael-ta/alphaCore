library(igraph)
library(dplyr)
library(fda.usc)
library(depthTools)
library(network)
library(car)
source("helper.R")


setwd("/mnt/alphaCore")

# set up graph if in list format
data<-read.csv(file="./data/networkcitation.txt", sep="\t", header=F)
#remove duplicates agfter april 2018, duplicates are due to a bug in java code.
data<-distinct(data)
#data=data[250:290,]#10 edges, 16 nodes
colnames(data)<-c("from","to","time","weight")
#create a directed multiple graph
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE) 
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE)%>%
    set_vertex_attr("idx", value = as.character(seq(1,vcount(tokenGr))))
#tokenGr2=network(data.frame(data[,1:2]),type="edgelist")#using "network" package, track node
#create weights
E(tokenGr)$weight<-data$weight


#Vertex count:
vcount(tokenGr)

#Igraph expects sequential ides, but our edge list ids start from 7M, igraph will create the missing ids
#We need to remove these artifical ids
#tokenGr<-delete_vertices((tokenGr), degree(tokenGr)==0)%>%#track node from original network, add attribute of node, i.e. origional index
#  set_vertex_attr("idx", value = as.character(seq(1,vcount(tokenGr))))

tokenGr<-delete_vertices((tokenGr), degree(tokenGr)==0)

vcount(tokenGr)
ecount(tokenGr)
idx_track=vertex_attr(tokenGr)
pdf("before-alphacore.pdf") 
#plot(tokenGr,vertex.size=2,vertex.label.cex=0.5,edge.arrow.size=0.3,rescale=FALSE,asp = 1.5)
# plot without labels for now
plot(tokenGr, vertex.size=2, vertex.label=NA, edge.arrow.size=0.1,rescale=TRUE)
dev.off() 

edge=as_edgelist(tokenGr, names = TRUE)
colnames(edge)<-c("Source","Target")
write.csv(edge, file = "edgelist.csv")

id=1:vcount(tokenGr)
node=cbind(id, V(tokenGr)$idx, 1:vcount(tokenGr) * 0)
colnames(node)<-c("Id","Label","group")
write.csv(node, file = "node.csv")

#Below are the alpha core steps

#when alpha>0, calculate depth value of each node and delete the nodes with depth value 
#>=alpha, for each alpha, run the iteration until no node is deleted under this alpha
#in each alpha, nodes are deleted;the update network is assigned with new depth value and comparred with the same alpha as before
#aCore returns alphaCoreMap which contains a rank of nodes from core to less core.
aCore<-function(tokenGr, alphaCoreMap,step=0.05){
  alphaCoreMap<-c();
  alpha<-1.0-step;
  m<-matrix();

  updated<-TRUE;
  level <- NULL
  min.weight <- NULL
  power.sample <- NULL

  while(alpha>=0){
    
    if(vcount(tokenGr)==0){
      message("Graph has no nodes left.")
      return(alphaCoreMap);
    }
    
    depthInputData<-data.frame();
    for(v in V(tokenGr)){#for each vertix in graph
      inDegree<-length(unlist(adjacent_vertices(tokenGr, v, mode = "in")))#count in degree
      inWeight<-sum(as.numeric(incident(tokenGr, v, mode = "in")$weight))#count weight of income edges
      #message("Node neighbor:",inDegree,", weight:",inWeight)
      newRow<-c(v,inDegree,inWeight)# node index, in degree, in weight
      depthInputData<-rbind(depthInputData,newRow)
    }

    # weight normalization
    #idx <- NULL
    #for (i in 1:max(depthInputData[,2])) {
    #    idx = which(depthInputData[,2] == i)
    #    if (is.null(min.weight)) {
    #        min.weight = min(depthInputData[idx,3])
    #    }
    #    if (length(idx) == 0) {
    #        next
    #    }
    #    depthInputData[idx,3] = depthInputData[idx, 3] - i * min.weight
    #}

    # transform data to be closer to normal
    #depthInputData[, 3] = depthInputData[, 3] ** .05
    #depthInputData[, 3] = depthInputData[, 3] ** .5
    if (is.null(power.sample)) {
        bcmodel <- boxCox(depthInputData[,3]~1, family="yjPower", plotit=F)
        power.sample = bcmodel$x[ which.max(bcmodel$y) ]
    }
    depthInputData[,3] = yjPower(depthInputData[,3], power.sample)


    #This depth function should be implemented
    #depthValue<-runif(nrow(depthInputData),0,1)
    colnames(depthInputData)<-c("node","inDegree","inWeight")
    write.csv(depthInputData, file = "depthInputData.csv")
    #calculate depth value of each node w.r.t. all other nodes
    #5 multivariate depth (Tukey(half space) depth, simplicial depth, Mahalanobis depth, random projeciton depth and likelihood depth)
    #2 functional depth will be used
    #Half space
    #temp=mdepth.HS(depthInputData[,2:3])
    #Simplicial depth
    #temp=mdepth.SD(depthInputData[,2:3])
    #MhD
    # calculate the mahalanobis distance w.r.t. the origini using indegree and inweight
    # keep sample variance consistent
    temp=1/(1 + mahalanobis.origin(depthInputData[,2:3], cov(depthInputData[,2:3])))
    #Random projection
    #temp=mdepth.RP(depthInputData[,2:3])
    #Likelihood depth
    #temp=mdepth.LD(depthInputData[,2:3])
    #MBD
    #temp=MBD(depthInputData[,2:3],plotting=TRUE)
    #ED
    
    #for functional depth
    #depthValue=as.matrix(temp$MBD)
    depthValue = temp
    depthInputData<-cbind(depthInputData, depthValue)
    colnames(depthInputData)<-c("node","inDegree","inWeight","depth")#last col: depth value
    rownames(depthInputData)<-NULL
    updated=FALSE;
    if (is.null(level)) {
        level <- sort(depthValue, decreasing=T)[vcount(tokenGr) * 0.2]
    }

    message("Alphacore is running for alpha:", alpha)
    nodesToBeDeleted <- which(depthValue >= level)
   
    if (length(nodesToBeDeleted) > 0) {
      updated=TRUE;
      message("Mean weight of deleted nodes: ", mean(depthInputData[nodesToBeDeleted,3]))
      alphaCoreMap <- rbind(alphaCoreMap, 
                            cbind(vertex_attr(tokenGr)$idx[nodesToBeDeleted], matrix(alpha, 
                                  length(nodesToBeDeleted))))
    }
    
    #for( v in V(tokenGr)){
      
      #depth comparison
    # if(depthInputData[depthInputData$node==v,]$depth>=alpha){
        #record alpha value for the node
        #alphaCoreMap<-rbind(alphaCoreMap,c(v,alpha))
    #    alphaCoreMap<-rbind(alphaCoreMap,c(vertex_attr(tokenGr)$idx[v],alpha))#record node index in origional network
    #    nodesToBeDeleted<-c(v, nodesToBeDeleted)
    #   updated=TRUE;
    #  }
    #}
    if(updated==FALSE){
      message("Nothing was removed for alpha: ", alpha);
      alpha = alpha-step;
      level = NULL;
    }
    else {
      message(length(nodesToBeDeleted), " nodes are deleted for alpha:",alpha,". Graph has ",vcount(tokenGr)," nodes.")
      tokenGr = delete_vertices(tokenGr, nodesToBeDeleted);
      level = 1;
    }
  }
    
  message(alpha," leaves with ", vcount(tokenGr))

  if (vcount(tokenGr) > 0) {
    for (v in V(tokenGr)) {
        alphaCoreMap<-rbind(alphaCoreMap,c(vertex_attr(tokenGr)$idx[v], 0))#record node index in origional network
    }
  } 

  return(alphaCoreMap)
}

#Store initial neighbor and weights as a record.
initialNodeFeatures<-data.frame()
for(v in V(tokenGr)){
  inDegree<-length(unlist(adjacent_vertices(tokenGr, v, mode = "in")))#count in degree
  inWeight<-sum(incident(tokenGr, v, mode = "in")$weight)
  #message("Node neighbor:",inDegree,", weight:",inWeight)
  newRow<-c(v,inDegree,inWeight)
  initialNodeFeatures<-rbind(initialNodeFeatures,newRow)
} 
colnames(initialNodeFeatures)<-c("node","inDegree","inWeight")

#run alpha core
alphaCoreMap<-aCore(tokenGr)

#plot network, with color of node indicates order of being deleted from the network


#correlation of core value and indegree and inweight
alphaCoreMap=cbind(1:vcount(tokenGr),alphaCoreMap)
colnames(alphaCoreMap)<-c("rank","node","alpha")

# buggy
vcolor <- c()
count <- 0
#print(sort(unique(as.numeric(alphaCoreMap[,3]))))
pdf("degree_vs_weights.pdf")
plot(initialNodeFeatures[,2:3], pch=".")

first <- TRUE
level <- 0
total <- length(unique(as.numeric(alphaCoreMap[,3])))
for (alpha in sort(unique(as.numeric(alphaCoreMap[,3])))) {
   level = level + 1
   idx  <- which(as.numeric(alphaCoreMap[,3]) == alpha)
   if (first) {
        print( initialNodeFeatures[ alphaCoreMap[idx,2], ] )
        first <- FALSE
   }
   count = count + length(idx)
   color <- rgb(1 - (count/vcount(tokenGr) * (level/total)), 0, count/vcount(tokenGr) * (level/total) )
   points(initialNodeFeatures[ alphaCoreMap[idx,2]  ,2:3], col=color, pch=".")
   vcolor[ alphaCoreMap[idx,2] ] = color
}
dev.off()

#for (i in 1:vcount(tokenGr)) {
#    alpha = as.numeric(alphaCoreMap[i, 3])
#    idx = as.numeric(alphaCoreMap[i, 2]) 
#    vcolor[idx] = rgb(1 - alpha, 0, alpha)
#} 

vertex_attr(tokenGr)$color = vcolor


pdf("after-alphacore.pdf") 
plot(tokenGr,vertex.size=2,vertex.label=NA,edge.arrow.size=0.1,rescale=T)
dev.off()

write.csv(alphaCoreMap, "results.csv")

#sort matrix according to node index
alphaCoreMap=alphaCoreMap[order(alphaCoreMap[,2]),]
temp=cbind(initialNodeFeatures,alphaCoreMap[,1],alphaCoreMap[,3])
colnames(temp)<-c("node","inDegree","inWeight","rank","alpha")
cor_indegree_rank=corr(temp[,2],temp[,4])
cor_inweight_rank=corr(temp[,3],temp[,4])
cor_indegree_alpha=corr(temp[,2],temp[,5])
cor_inweight_alpha=corr(temp[,3],temp[,5])




