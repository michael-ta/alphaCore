  library(igraph)
  library(qgraph)
  library(dplyr)
  library(fda.usc)
  library(depthTools)
  library(network)
  library(car)
  
  args <- commandArgs(trailing=F)
  script.path <- sub("--file=", "", args[grep("--file=", args)])
  if (!(identical(script.path, character(0)))) {
    script.basename <-  dirname(script.path)
  } else {
    # for debug / rstudio purposes
    script.basename <- "/mnt/alphaCore"
  }
  setwd(script.basename) 
  source("helper.R")
  
  args <- commandArgs(trailing=T)
  # debug parameters
  debug.network <- "airport-US2010"
  debug.stepsize <- 0.05
  data.network <- if (!is.na(args[1])) args[1] else debug.network
  step.size <- if (!is.na(args[2])) as.numeric(args[2]) else debug.stepsize
  
  data.network.fn <- paste("./data/network", data.network, "txt", sep=".")
  data.labels.fn <- paste("./data/label", data.network, "txt", sep=".")
  data.labels.fn <- if (!(file.exists(data.labels.fn))) NA else data.labels.fn
 
  data<-read.csv(file=data.network.fn, sep=" ", header=F)
  data[,3] <- (data[,3] - min(data[,3])) / (max(data[,3]) - min(data[,3]))
  
  if (is.na(data.labels.fn)) {
      data.labels <- NULL
  } else {
      data.labels<-read.csv(file=data.labels.fn, header=F)
      data.labels<-as.character(data.labels$V1)
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
  #write.csv(edge, file = "edgelist.csv")
  
  id=1:vcount(tokenGr)
  node=cbind(id, V(tokenGr)$idx, 1:vcount(tokenGr) * 0)
  colnames(node)<-c("Id","Label","group")
  #write.csv(node, file = "node.csv")
  
  # Below are the alpha core steps
  #
  # when alpha>0, for mahalanobis depth 
  #   calculate depth value of each node w.r.t the origin and delete 20% 
  #   (parameter) of nodes with the highest depth value. If nodes were deleted 
  #   in the current iteration, re-run and delete nodes with 1 depth i.e. nodes
  #   without any ingoing edges. Returns alphaCoreMap which contains a rank of 
  #   nodes from core to less core.
  aCore<-function(tokenGr, edgelist, depthInputData, alphaCoreMap,step=0.01){
    alphaCoreMap<-c();
    alpha<-1
  
    updated<-TRUE;
    level <- NULL
    min.weight <- NULL
    min.depth <- NULL
    max.degree <- max(depthInputData[,2])
    depthInputData <- cbind(depthInputData, 0, 0, 1) 
    
    while(alpha > 0){
      if(vcount(tokenGr)==0){
        message("Graph has no nodes left.")
        return(alphaCoreMap);
      }
      
      # normalized the degree - we've already normalize the weights
	    depthInputData[,5] = depthInputData[,3]
	    depthInputData[,4] = depthInputData[,2] / max.degree
  
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
                      depthInputData[,4:5], 
                      cov(depthInputData[,4:5])))
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
      colnames(depthInputData)<-c("node",
                                  "inDegree",
                                  "inWeight", 
                                  "tr_inDegree", 
                                  "tr_inWeight", 
                                  "depth")
      rownames(depthInputData)<-NULL
      updated=FALSE;
  
      # get the depth value to remove if a new alpha level is reached 
      # (otherwise we set to 1 from prior iteration after removing nodes 
      # at prevoius alpha level)
      if (is.null(level)) {
          level <- alpha 
      }
  
      message("Alphacore is running for alpha:", alpha)
      nodesToBeDeleted <- which(depthInputData[,6] >= level)
     
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
          eidx <- which(edgelist[,1] == v)
          if (!(is.null(nrow(depthInputData)))) {
              didx <- which(depthInputData[,1] == v)
              for (anode in edgelist[eidx,2]) {
                  aidx <- which(depthInputData[,1] == anode)
                  depthInputData[aidx,2] = depthInputData[aidx, 2] - 1
                  widx <- which(edgelist[eidx,2] == anode)
                  depthInputData[aidx,3] = depthInputData[aidx, 3] - 
                                             sum(edgelist[eidx,3][widx])
              }
              depthInputData = depthInputData[-didx,]
          } else {
              depthInputData = NULL
              alpha = 0
          }
  
          deleted_tmp = c(deleted_tmp, 
                          which(as.integer(vertex_attr(tokenGr)$idx) == v))
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
  alphaCoreMap<-aCore(tokenGr, data, initialNodeFeatures, step=step.size)
  
  # plot network, where the color of node indicates the order of being deleted 
  # from the network, from furthest to closest blue -> purple -> red -> green
  alphaCoreMap=cbind(1:vcount(tokenGr),alphaCoreMap)
  colnames(alphaCoreMap)<-c("rank","node","alpha")
  rownames(alphaCoreMap)<-c()
  alphaCoreMap.labels <- c()

  
  write.csv(cbind(alphaCoreMap, alphaCoreMap.labels), 
            paste("results", data.network, 
                  format(step.size, scientific=T), "csv", sep="."))
  
  top.nodes.idx <- which(initialNodeFeatures[,1] %in% 
                         alphaCoreMap[which(as.numeric(alphaCoreMap[,3]) == 
                                            min(alphaCoreMap[,3])), 2])
  top.weights <- sum(initialNodeFeatures[top.nodes.idx,3]) / 
                    sum(initialNodeFeatures[,3])
  top.edges <- sum(initialNodeFeatures[top.nodes.idx, 2]) /
                    sum(initialNodeFeatures[,2])

  
  message("\n\n----- Alphacore Identifies -----")
  message("Top ", length(top.nodes.idx), " nodes account for ", 
          sprintf("%.3f", top.weights * 100), "% of network weights")

  message("Top ", length(top.nodes.idx), " nodes account for ", 
          sprintf("%.3f", top.edges *100), "% of network edges")
  message("---------------------------------")

  #sum(initialNodeFeatures[which(initialNodeFeatures[,1] %in% top.nodes.idx),3]) / sum(initialNodeFeatures[,3])

  #cor_indegree_rank=cor(temp[,2],temp[,4])
  #cor_inweight_rank=cor(temp[,3],temp[,4])
  #cor_indegree_alpha=cor(temp[,2],temp[,5])
  #cor_inweight_alpha=cor(temp[,3],temp[,5])
