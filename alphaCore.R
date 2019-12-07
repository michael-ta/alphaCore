library(igraph)
library(qgraph)
library(dplyr)
library(fda.usc)
library(depthTools)
library(network)
library(car)
source("helper.R")


setwd("/mnt/alphaCore")
data.idx <- 2

data.fn <- c("./data/networkcitation.txt",
             "./data/US_airport_2010.txt",
             "./data/4932.protein.links.v11.0.small.txt")

data.labels.fn <- c("./data/authorList.txt",
                    "./data/USairport_2010_codes.txt",
                    NULL)

data<-read.csv(file=data.fn[data.idx], sep=" ", header=F)

if (is.null(data.labels.fn)) {
    data.labels <- NULL
} else {
    data.labels<-read.csv(file=data.labels.fn[data.idx], header=F)
    data.labels<-as.character(data.labels$V1)
}

# remove duplicates agfter april 2018, duplicates are due to a bug in java 
# code.
data<-distinct(data)
colnames(data)<-c("from", "to", "weight")
# create a directed graph and store the initial node index in the vertex 
# attribute idx
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE) 
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE)%>%
    set_vertex_attr("idx", value = as.character(seq(1,vcount(tokenGr))))

# add weights to each edge in the graph
E(tokenGr)$weight<-data$weight

# delete vertices that iGraph fills in, or any vertices that have 0 degrees 
# i.e. not connected to the rest of the network
tokenGr<-delete_vertices((tokenGr), degree(tokenGr)==0)

# simplify graph to remove self loops
tokenGr<-simplify(tokenGr)

vcount(tokenGr)
ecount(tokenGr)

pdf("before-alphacore.pdf") 
plot(tokenGr, vertex.size=2, vertex.label=NA, edge.arrow.size=0.1,rescale=TRUE)
dev.off() 

edge=as_edgelist(tokenGr, names = TRUE)
colnames(edge)<-c("Source","Target")
write.csv(edge, file = "edgelist.csv")

id=1:vcount(tokenGr)
node=cbind(id, V(tokenGr)$idx, 1:vcount(tokenGr) * 0)
colnames(node)<-c("Id","Label","group")
write.csv(node, file = "node.csv")


# function to get edge weights
getEdgeWeights<-function(inputGr) {
    depthInputData<-data.frame();
    # for each vertex in graph
    for(v in V(inputGr)) {
        # count in degree
        inDegree<-length(unlist(adjacent_vertices(inputGr, v, mode="in")))
        # count weight of incoming edges
        inWeight<-sum(as.numeric(incident(inputGr, v, mode="in")$weight))
        # node index, in degree, in weight
        newRow<-c(v, inDegree, inWeight)
        depthInputData<-rbind(depthInputData, newRow)
    }
    depthInputData
}

# function to color nodes given the graph and alphaCoreMap output from 
# alphaCore
getVertexColors<-function(inputGr, aCM) {
    vcolor <- c()
    total <- length(unique(as.numeric(aCM[,3])))
    level <- 0
    count <- 0
    first <- TRUE
    for (alpha in sort(unique(as.numeric(aCM[,3])), decreasing=T)) {
        level = level + 1
        idx  <- which(as.numeric(aCM[,3]) == alpha)
        count = count + length(idx)
        color <- rgb(1 - (count/vcount(inputGr) * (level/total)), 
                     0,
                     count/vcount(inputGr) * (level/total) )
        tidx <- which(as.numeric(vertex_attr(inputGr, "idx")) %in% 
                      as.numeric(aCM[idx,2]))

        if (first) {
            color <- rgb(0, 1, 0)
            first <- FALSE
        }
        vcolor[ tidx ] = color
    }
    vcolor
}

# Below are the alpha core steps
#
# when alpha>0, for mahalanobis depth 
#   calculate depth value of each node w.r.t the origin and delete 20% 
#   (parameter) of nodes with the highest depth value. If nodes were deleted 
#   in the current iteration, re-run and delete nodes with 1 depth i.e. nodes
#   without any ingoing edges. Returns alphaCoreMap which contains a rank of 
#   nodes from core to less core.
aCore<-function(tokenGr, alphaCoreMap,step=0.01){
  alphaCoreMap<-c();
  #alpha<-1.0-step;
  alpha<-1
  m<-matrix();

  updated<-TRUE;
  level <- NULL
  min.weight <- NULL
  power.sample <- NULL
  step.est <- NULL

  #while(alpha>=0){
  while(TRUE){
    
    if(vcount(tokenGr)==0){
      message("Graph has no nodes left.")
      return(alphaCoreMap);
    }
    
    depthInputData<-getEdgeWeights(tokenGr);
    g.density <- density(depthInputData[,2])
    rmCount = max(which(g.density[[1]] <= 1))
    plot(g.density)
    points(g.density[[1]][rmCount], g.density[[2]][rmCount], col="red")

    # perform a power transformation on the variables
    if (is.null(power.sample)) {
        bcmodel <- boxCox(depthInputData[,3]~1, family="yjPower", plotit=F)
        power.sample.weight = bcmodel$x[ which.max(bcmodel$y) ]
        bcmodel <- boxCox(depthInputData[,2]~1, family="yjPower", plotit=F)
        power.sample.degree = bcmodel$x[ which.max(bcmodel$y) ]
        #message("Using power.sample.degree: ", power.sample.degree)
        #message("Using power.sample.weight: ", power.sample.weight)
    }
    depthInputData[,3] = yjPower(depthInputData[,3], power.sample.weight)
    depthInputData[,2] = yjPower(depthInputData[,2], power.sample.degree)

    if (is.null(step.est)) {
        step.est = sum(depthInputData[,2] == min(depthInputData[,2])) / vcount(tokenGr) 
        message("Using step size: ", step.est)
    }

    #This depth function should be implemented
    #calculate depth value of each node w.r.t. all other nodes
    #5 multivariate depth (Tukey(half space) depth, simplicial depth, Mahalanobis depth, random projeciton depth and likelihood depth)
    #2 functional depth will be used
    #Half space
    #temp=mdepth.HS(depthInputData[,2:3])
    #Simplicial depth
    #temp=mdepth.SD(depthInputData[,2:3])
    #MhD
    # calculate the mahalanobis depth w.r.t. the origin using indegree and 
    # inweight keep sample variance consistent
    temp=1/(1 + mahalanobis.origin(
                    depthInputData[,2:3], 
                    cov(depthInputData[,2:3])))
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
    # add column for depth to dataframe
    depthInputData<-cbind(depthInputData, depthValue)
    colnames(depthInputData)<-c("node","inDegree","inWeight","depth")
    rownames(depthInputData)<-NULL
    updated=FALSE;

    # get the depth value associated with the 20 percentile of all nodes in 
    # the current graph
    if (is.null(level)) {
        level <- sort(depthValue, decreasing=T)[ ceiling(vcount(tokenGr) * step.est)]
        level <- sort(depthValue, decreasing=T)[ ceiling(vcount(tokenGr) * .10)]
        #level <- sort(depthValue, decreasing=T)[rmCount]
    }

    message("Alphacore is running for alpha:", alpha)
    nodesToBeDeleted <- which(depthInputData[,4] >= level)
   
    if (length(nodesToBeDeleted) > 0) {
      updated=TRUE;
      # record alpha level for deleted nodes prior to removing them
      alphaCoreMap <- rbind(alphaCoreMap, 
                            cbind(vertex_attr(tokenGr)$idx[nodesToBeDeleted],
                                  matrix(alpha, length(nodesToBeDeleted))))
    }
    
   if(updated==FALSE){
      message("Nothing was removed for alpha: ", alpha);
      # fix rounding errors for decimal values after many iterations
      #alpha = signif(alpha-step, 2);
      alpha = alpha + 1;
      if (is.na(level) || level != 1) {
          break
      }
      level = NULL;
      power.sample.weight = 0;
      power.sample.degree = 0;
    }
    else {
      message(length(nodesToBeDeleted), 
              " nodes are deleted for alpha:", alpha, 
              ". Graph has ",vcount(tokenGr)," nodes.")
      tokenGr = delete_vertices(tokenGr, nodesToBeDeleted);
      level = 1;
    }
  }
    
  message(alpha," leaves with ", vcount(tokenGr))

  if (vcount(tokenGr) > 0) {
    for (v in V(tokenGr)) {
        alphaCoreMap<-rbind(alphaCoreMap,c(vertex_attr(tokenGr)$idx[v], alpha))
        #record node index in origional network
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

# plot network, where the color of node indicates the order of being deleted 
# from the network, from furthest to closest blue -> purple -> red -> green
alphaCoreMap=cbind(1:vcount(tokenGr),alphaCoreMap)
colnames(alphaCoreMap)<-c("rank","node","alpha")

vcolor <- c()
count <- 0
pdf("degree_vs_weights.pdf")
plot(initialNodeFeatures[,2:3], pch=".")

first <- TRUE
level <- 0
total <- length(unique(as.numeric(alphaCoreMap[,3])))
alphaCoreMap.labels <- c()

for (alpha in sort(unique(as.numeric(alphaCoreMap[,3])), decreasing=T)) {
   level = level + 1
   idx  <- which(as.numeric(alphaCoreMap[,3]) == alpha)

   count = count + length(idx)
   color <- rgb(1 - (count/vcount(tokenGr) * (level/total)), 
                0,
                count/vcount(tokenGr) * (level/total) )
   # incorrect previously used tokenGr vertex attr
   tidx <- which(as.numeric(vertex_attr(tokenGr, "idx")) %in% as.numeric(alphaCoreMap[idx,2]))

   if (first) {
        print( cbind(alphaCoreMap[idx, 2], initialNodeFeatures[ tidx, 2:3 ]) )
        color <- rgb(0, 1, 0)
        first <- FALSE
   }

   #alphaCoreMap.labels = c(as.matrix(
   #                        data.labels)[as.numeric(alphaCoreMap[idx, 2])],
   #                        alphaCoreMap.labels)

   points(initialNodeFeatures[ tidx  ,2:3], col=color, pch=".")
   vcolor[ tidx ] = color
}
dev.off()

vertex_attr(tokenGr)$color = vcolor
vertex_attr(tokenGr)$label.cex = c(rep(0.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.dist = c(rep(.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.degree = c(rep(-pi/6, vcount(tokenGr)))
vertex_attr(tokenGr)$width = c(rep(0.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.color = c(rep("black", vcount(tokenGr)))

tokenGr <- tokenGr %>% set_edge_attr("color", value=rgb(0.7, 0.7, 0.7, 0.25))

didx <- which(as.numeric(alphaCoreMap[,3]) >= max(as.numeric(alphaCoreMap[,3])) - 10)
dnodes <- as.numeric( alphaCoreMap[didx, 2] )
# get original node ids of dnodes for subgraph creation
sidx <- which(as.numeric(vertex_attr(tokenGr)$idx) %in% dnodes)
sGr <- induced_subgraph(tokenGr, sidx)



pdf("remaining-alphacore.pdf", width=12, height=12)
e <- get.edgelist(sGr)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(sGr))
if (is.null(data.labels)) {
    labels=NA
} else {
    labels=data.labels[sort(dnodes)]
}

vertex_attr(sGr)$color = getVertexColors(sGr, alphaCoreMap[didx,])

plot(sGr, 
     layout=l, 
     vertex.size=2,
     label.dist=.5,
     label.cex=0.05, 
     vertex.label=labels,
     edge.arrow.size=0.1,
     width=0.25,
     rescale=T, asp=0)
dev.off()



pdf("after-alphacore.pdf", width=24, height=24)
e <- get.edgelist(tokenGr)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(tokenGr)) 
plot(tokenGr, 
     layout=l,
     vertex.size=2,
     vertex.label=NA,
     edge.arrow.size=0.1,
     rescale=T, asp=0)
dev.off()

write.csv(cbind(alphaCoreMap, alphaCoreMap.labels), "results.csv")

#sort matrix according to node index
alphaCoreMap=alphaCoreMap[order(alphaCoreMap[,2]),]
temp=cbind(initialNodeFeatures,alphaCoreMap[,1],alphaCoreMap[,3])
colnames(temp)<-c("node","inDegree","inWeight","rank","alpha")

# correlation of core value and indegree and inweight
# find correlations (?)
cor_indegree_rank=cor(temp[,2],temp[,4])
cor_inweight_rank=cor(temp[,3],temp[,4])
cor_indegree_alpha=cor(temp[,2],temp[,5])
cor_inweight_alpha=cor(temp[,3],temp[,5])
