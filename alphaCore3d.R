  library(igraph)
  library(qgraph)
  library(dplyr)
  library(fda.usc)
  library(depthTools)
  library(network)
  library(car)
  
# Version V.2.2
# incorporation of generalized MhD for nominal data
  
  setwd("/mnt/alphaCore")
  source("helper.R")
  
  data.idx <- 3
  step.size <- 0.000005
  
  data.fn <- c("./data/network.test-network.txt",
               "./data/network.metal-trade.txt",
               "./data/network.BAT.txt")
  
  data.labels.fn <- c("./data/label.test-network.txt",
                      "./data/label.metal-trade.txt",
                      NA)

  data.features.fn <- c("./data/features.test-network.txt",
                        "./data/features.metal-trade.txt",
                        "./data/features.BAT.txt")
  
  data<-read.csv(file=data.fn[data.idx], sep=" ", header=F)
  # normalize the input weights
  data[,3] <- (data[,3] - min(data[,3])) / (max(data[,3]) - min(data[,3]))
  
  if (is.na(data.labels.fn[data.idx])) {
      data.labels <- NULL
  } else {
      data.labels<-read.csv(file=data.labels.fn[data.idx], header=F)
      data.labels<-as.character(data.labels$V1)
  }

  if (is.na(data.features.fn[data.idx])) {
      data.features <- NULL
  } else {
      # format of features should be 
      #     node_id, feature_1, feature_2, ..., feature_p
      data.features <- read.csv(file=data.features.fn[data.idx], 
                                header=T, sep=" ")  
  }
  
  # remove duplicates agfter april 2018, duplicates are due to a bug in java 
  # code.
  data <- distinct(data)
  # remove self-loops
  duplicates <- which(data[,1] == data[,2])
  if (!(length(duplicates) == 0)) {
    data = data[-duplicates,]
  }
  colnames(data)<-c("from", "to", "weight")
  # create a directed graph and store the initial node index in the vertex 
  # attribute idx
  tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE) 
  #tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE)%>%
  tokenGr<-tokenGr%>%
      set_vertex_attr("idx", value = as.character(seq(1,vcount(tokenGr))))
  
  # add weights to each edge in the graph
  E(tokenGr)$weight<-data$weight
  
  # delete vertices that iGraph fills in, or any vertices that have 0 degrees 
  # i.e. not connected to the rest of the network
  tokenGr<-delete_vertices((tokenGr), degree(tokenGr)==0)
  
  # simplify graph to remove self loops
  #tokenGr<-simplify(tokenGr)
  
  vcount(tokenGr)
  ecount(tokenGr)
  
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
      for (alpha in sort(unique(as.numeric(aCM[,3])), decreasing=F)) {
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
 

  transform.BC <- function(x) {
    bcmodel <- boxCox(x~1, family="yjPower", plotit=F)
    power.sample <- bcmodel$x[which.max(bcmodel$y) ]
    result = yjPower(x, power.sample)
    return(result)
  }
 
  # Below are the alpha core steps
  #
  # when alpha>0, for mahalanobis depth 
  #   calculate depth value of each node w.r.t the origin and delete 20% 
  #   (parameter) of nodes with the highest depth value. If nodes were deleted 
  #   in the current iteration, re-run and delete nodes with 1 depth i.e. nodes
  #   without any ingoing edges. Returns alphaCoreMap which contains a rank of 
  #   nodes from core to less core.
  # params
  #   features - dataframe with node features (id, v1, ..., vp), if feature is
  #              nominal it needs to be encoded with as.factor - otherwise it
  #              is assumed to be continuous and we will join the featue with 
  #              the degree / weight matrix

  aCore<-function(tokenGr, edgelist, features,  depthInputData, alphaCoreMap,step=0.01){
    alphaCoreMap<-c();
    alpha<-1
  
    updated<-TRUE;
    level <- NULL
    min.weight <- NULL
    power.sample <- NULL
    min.depth <- NULL
    max.degree <- max(depthInputData[,2])
    depthInputData <- cbind(depthInputData, 0, 0, 1)
 
    # index of factor features
    faidx <- grep("factor", colnames(features))
    # index of continuous features
    ctidx <- grep("count", colnames(features))  
  
    while(alpha > 0){
      if(vcount(tokenGr)==0){
        message("Graph has no nodes left.")
        return(alphaCoreMap);
      }
      
      # normalized the degree - we've only normalize the weights
      depthInputData[,5] = depthInputData[,3]
      depthInputData[,4] = depthInputData[,2] / max.degree

      # calculate MhDepth for the factor features (currently unused)
      #mhd.factor <- mahalanobis.factor(features)
      # set this to node id
      #mhd.factor[,1] = features[,1]     
  
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
      depthInputData_tmp <- merge(depthInputData, features[,c(1, ctidx)], by="node")

      # rescale by factors (currently unused)
      #temp=1/(1 + mahalanobis.origin(
      #                depthInputData_tmp[,c(4,5,7:(6+length(ctidx)))], 
      #                cov(depthInputData_tmp[,c(4,5,7:(6+length(ctidx)))])) * mhd.factor[,2] + mhd.factor[,1]  )

      temp=1/(1 + mahalanobis.origin(
                      depthInputData_tmp[,c(4,5,7:(6+length(ctidx)))], 
                      cov(depthInputData_tmp[,c(4,5,7:(6+length(ctidx)))])) )

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
      depthInputData[,6]<- depthValue
      colnames(depthInputData)<-c("node","inDegree","inWeight", "tr_inDegree", "tr_inWeight", "depth")
      rownames(depthInputData)<-NULL
      updated=FALSE;
  
      # get the depth value associated with the 20 percentile of all nodes in 
      # the current graph
      if (is.null(level)) {
          level <- alpha 
      }
  
      message("Alphacore is running for alpha:", alpha)
      # TODO:
      # for 3D alphaCore we need to additionally remove nodes that do not have
      # any incoming edges even though the depth value != 0. This is because
      # the node features are non-zero

      nodesToBeDeleted <- which(depthInputData[,2] == 0)
      nodesToBeDeleted <- union(nodesToBeDeleted, which(depthInputData[,6] >= level))

     
      if (length(nodesToBeDeleted) > 0) {
        updated=TRUE;
        # record alpha level for deleted nodes prior to removing them
        nodesToBeDeleted <- depthInputData[nodesToBeDeleted, 1]
        alphaCoreMap <- rbind(alphaCoreMap, 
                              cbind(nodesToBeDeleted,
                                    matrix(alpha, length(nodesToBeDeleted))))
      }
      
     if(updated==FALSE){
        message("Nothing was removed for alpha: ", alpha);
        # fix rounding errors for decimal values after many iterations
        delta = (alpha - sort(depthInputData[which(depthInputData[,6] < alpha),6], decreasing=T)[1]) / step
        prev_alpha_tmp = alpha
        alpha = signif(alpha - (step * ceiling(delta)), 15)
        if (is.na(alpha)) {
          alpha = 0
          # alphaCore reiterate
          #delta = (prev_alpha_tmp - sort(depthInputData[,6], decreasing=T)[1]) / step
          #alpha = signif(prev_alpha_tmp - (step * ceiling(delta)), 15)
          #if (is.na(alpha)) {
          #  alpha = 0
          #}
        }
       level = NULL;
        power.sample.weight = 0;
        power.sample.degree = 0;
      }
      else {
        message(length(nodesToBeDeleted), 
                " nodes are deleted for alpha:", alpha, 
                ". Graph has ",vcount(tokenGr)," nodes.")
        # convert nodesToBeDeleted to match vertex id
        # recalculate depthInputData degree and weights for deleted nodes
        deleted_tmp = c()
        for (v in nodesToBeDeleted) {
          # affected edges if we remove this node
          eidx <- which(edgelist[,1] == v)
          if (!(is.null(nrow(depthInputData)))) {
              # index of inweight, indegree
              didx <- which(depthInputData[,1] == v)
              # iterate through affected nodes if we make this deletion
              for (anode in edgelist[eidx,2]) {
                  aidx <- which(depthInputData[,1] == anode)
                  # update degree
                  depthInputData[aidx,2] = depthInputData[aidx, 2] - 1
                  widx <- which(edgelist[eidx,2] == anode)
                  depthInputData[aidx,3] = depthInputData[aidx, 3] - edgelist[eidx,3][widx]
              }
              # remove node from node list
              depthInputData = depthInputData[-didx,]
              # remove row from feature
              features = features[-which(features[,1] == v),]
          } else {
              depthInputData = NULL
              alpha = 0
          }
  
          deleted_tmp = c(deleted_tmp, which(as.integer(vertex_attr(tokenGr)$idx) == v))
        }
        tokenGr = delete_vertices(tokenGr, deleted_tmp);
        if (vcount(tokenGr) == 1) {
          alpha = 0
        }
        level = 1;
      }
    }
      
    message(alpha," leaves with ", vcount(tokenGr))
  
    if (vcount(tokenGr) > 0) {
      for (v in V(tokenGr)) {
          alphaCoreMap<-rbind(alphaCoreMap,c(vertex_attr(tokenGr)$idx[v], 0))
          #record node index in origional network
      }
    } 
  
    return(alphaCoreMap)
  }
  
  #Store initial neighbor and weights as a record.
  initialNodeFeatures<-matrix(0, vcount(tokenGr), 3)
  
  # should be faster to iterate though the dataset
  counter = 1
  for(v in as.integer(vertex_attr(tokenGr)$idx)) {
      nidx <- which(data[,2] == v)
      inDegree <- length(nidx)
      inWeight <- sum(data[nidx, 3])
      initialNodeFeatures[counter,] = c(v,inDegree, inWeight) 
      counter = counter + 1
  }
  
  colnames(initialNodeFeatures)<-c("node","inDegree","inWeight")
  
  #run alpha core
  alphaCoreMap<-aCore(tokenGr, data, data.features, initialNodeFeatures, step=step.size)
  
  # plot network, where the color of node indicates the order of being deleted 
  # from the network, from furthest to closest blue -> purple -> red -> green
  alphaCoreMap=cbind(1:vcount(tokenGr),alphaCoreMap)
  colnames(alphaCoreMap)<-c("rank","node","alpha")
  rownames(alphaCoreMap)<-c()
  vcolor <- c()
  count <- 0
  pdf(paste("degree_vs_weights", step.size, "pdf", sep="."))
  plot(initialNodeFeatures[,2:3], pch=".")
  
  first <- TRUE
  level <- 0
  total <- length(unique(as.numeric(alphaCoreMap[,3])))
  alphaCoreMap.labels <- c()
  
  for (alpha in sort(unique(as.numeric(alphaCoreMap[,3])), decreasing=F)) {
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
  vertex_attr(tokenGr)$size = 1 + ((initialNodeFeatures[,2] / max(initialNodeFeatures[,2])) * 3)

  tokenGr <- tokenGr %>% set_edge_attr("color", value=rgb(0.7, 0.7, 0.7, 0.25))
  
  didx <- c()
  for (i in 1:10) {
    new_didx <- which( as.numeric(alphaCoreMap[,3]) == sort(as.numeric(unique(alphaCoreMap[,3])))[i] )
    if (length(new_didx) + length(didx) < 3500) {
      didx = c(didx, new_didx)
    } else {
      didx = c(didx, new_didx[1:(3500 - length(didx))])
      break
    }
  }

  dnodes <- as.numeric( alphaCoreMap[didx, 2] )
  # get original node ids of dnodes for subgraph creation
  sidx <- which(as.numeric(vertex_attr(tokenGr)$idx) %in% dnodes)
  sGr <- induced_subgraph(tokenGr, sidx)
  
  
  pdf(paste("remaining-alphacore", step.size, "pdf", sep="."), width=12, height=12)
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
       label.dist=.5,
       label.cex=0.05, 
       vertex.label=labels,
       edge.arrow.size=0.1,
       width=0.25,
       rescale=T, asp=0)
  dev.off()
  
  write.csv(cbind(alphaCoreMap, alphaCoreMap.labels), paste("results", step.size, "csv", sep="."))
  
  #sort matrix according to node index
  alphaCoreMap=alphaCoreMap[order(alphaCoreMap[,2]),]
  temp=cbind(initialNodeFeatures,alphaCoreMap[,1],alphaCoreMap[,3])
  colnames(temp)<-c("node","inDegree","inWeight","rank","alpha")
  

  alevels <- as.data.frame(alphaCoreMap) %>% group_by(alpha) %>% summarise(no_rows = length(alpha))
  oidx <- order(as.numeric(alevels$alpha), decreasing=T)
  alevels.core <- cbind(cumsum(alevels$no_rows[oidx]) / sum(alevels$no_rows),
                        alevels$alpha[oidx])
  plot(alevels.core[,1], alevels.core[,2], type="b", pch=5,
       xlab="percentage of nodes", ylab="alpha level")

  top.nodes.idx <- which(initialNodeFeatures[,1] %in% alphaCoreMap[which(as.numeric(alphaCoreMap[,3]) == min(alphaCoreMap[,3])), 2])
  sum(initialNodeFeatures[top.nodes.idx,3]) / sum(initialNodeFeatures[,3])

  #sum(initialNodeFeatures[which(initialNodeFeatures[,1] %in% top.nodes.idx),3]) / sum(initialNodeFeatures[,3])

  # correlation of core value and indegree and inweight
  # find correlations (?)
  #cor_indegree_rank=cor(temp[,2],temp[,4])
  #cor_inweight_rank=cor(temp[,3],temp[,4])
  #cor_indegree_alpha=cor(temp[,2],temp[,5])
  #cor_inweight_alpha=cor(temp[,3],temp[,5])
