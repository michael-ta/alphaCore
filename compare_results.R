setwd("D:/repos/alphaCore/results/vary-stepsize")
filelist <- list.files(pattern="results.*.csv")

datalist = lapply(filelist, function(x) read.table(x, header=T, sep=","))

# these ids corresond to c("ATL", "LAX", "DFW", "SFO", "LAS", "PHX", "CLT", "MIA", "MSP")
toplist <- c(114, 877, 391, 389, 1437, 875, 1255, 317, 1016, 1068)

count = c()
network.size <- c()

for (i in 1:length(datalist)) {
  idx <- which(datalist[[i]]$alpha == max(datalist[[i]]$alpha))
  l <- length(datalist[[i]][idx,1])
  count = c(count, sum(datalist[[i]][idx,]$node %in% toplist))
  network.size = c(network.size, l)
}

plot(seq(0.01, 0.99, 0.01), count/10, type="l", xlab="Step size (0.01 increments)", ylab="# of top 10 airports in final iteration")
points(seq(0.01, 0.99, 0.01), count / network.size ,type="l", col="red")


### Find Kendall's Tau for groups with shared ranks

assign.ranks <- function(d) {
  prev_alpha = NA
  r = 0
  ranks = matrix(0, nrow(d), 1)
  for (idx in 1:nrow(d)) {
    if (is.na(prev_alpha)) {
      r = 1
      ranks[idx] = r
      prev_alpha = d$alpha[idx]
    } else if (d$alpha[idx] == prev_alpha) {
      ranks[idx] = r
    } else {
      r = r + 1
      ranks[idx ] = r
      prev_alpha = d$alpha[idx]
    }
  }
  cbind(d, ranks)
}


# implementation of Kendall Tau
# https://support.sas.com/documentation/cdl/en/procstat/63104/HTML/default/viewer.htm#procstat_corr_sect015.htm
kendall.tau <- function(d) {
  N = nrow(d)
  nCD = matrix(0, N-1, 2)
  
  for (idx in 1:N-1) {
    nc_tmp <- sum(d$rank.y[idx:N] > d$rank.y[idx])
    nd_tmp <- sum(d$rank.y[idx:N] < d$rank.y[idx])
    nCD[idx,] = c(nc_tmp, nd_tmp)
  }
  rX <- d %>% group_by(rank.x) %>% summarise(c=length(rank.x))
  rY <- d %>% group_by(rank.y) %>% summarise(c=length(rank.y))
  
  denom <- (N * (N-1) /2 - sum(rX$c * (rX$c - 1) / 2)) * (N * (N-1) /2 - sum(rY$c * (rY$c - 1) / 2))
  diff(rev(colSums(nCD))) / sqrt(denom)
}

setwd("D:/repos/alphaCore/results/v.2-Results/")
# Read in data and assign ranks

alpha <- c( "5e-02", "5e-03", "5e-04", "5e-05", "5e-06", "5e-07", "5e-08", "5e-09", "5e-10")

for(d in alpha[1:6]) {
datalist.1 <- read.table("richClub.netscience.csv", header=T, sep=",")
datalist.1 <- datalist.1[,2:4]
datalist.2 <- read.table(paste("netscience.reiterate/results.", d, ".csv", sep=""), header=T, sep=",")
datalist.2 <- datalist.2[,2:4]

#for rich club
datalist.1 <- datalist.1[order(datalist.1$bin),]
colnames(datalist.1) <- c("id", "node", "alpha")

datalist.1 <- assign.ranks(datalist.1)
colnames(datalist.1) <- c("id", "node", "alpha", "rank")
datalist.2 <- assign.ranks(datalist.2)
colnames(datalist.2) <- c("id", "node", "alpha", "rank")

datalist <- merge(datalist.1, datalist.2, by="node")
datalist <- datalist[order(datalist$rank.x),]
datalist <- datalist[,c("id.x", "node", "alpha.x", "alpha.y", "rank.x", "rank.y")]

print(kendall.tau(datalist))
plot(1:nrow(datalist.2)/nrow(datalist.2),datalist.2$alpha, pch=20, 
     ylab="alpha / MhD", 
     xlab="Percentage of Nodes", 
     main=paste("alphaCore Map reiter. (airport ", d, " )", collapse=""))
}


