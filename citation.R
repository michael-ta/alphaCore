library(igraph)

paperCitAdj <- read.table("data/Citation/paperCitAdj.txt")
authorPaperBiAdj <- read.table("data/Citation/authorPaperBiadj.txt")
authorList <- read.table("data/Citation/authorList.txt")

authorCitations <- matrix(0, nrow=nrow(authorPaperBiAdj), ncol=nrow(authorPaperBiAdj))

for (i in 1:nrow(paperCitAdj)) {
    # get the authors of current paper
    citedauthors <- which(authorPaperBiAdj[, i] == 1)
    # get papers that have been cited by others
    pidx <- which(paperCitAdj[i, ] == 1)
    if (length(pidx) > 0) {
        for (p in 1:length(pidx)) {
            authors <- which(authorPaperBiAdj[, pidx[p]] == 1)
            for (c in 1:length(citedauthors)) {
                authorCitations[citedauthors[c],] = authorCitations[citedauthors[c],] + t(authorPaperBiAdj[, pidx[p]])
            }
        }
    }
    if (i %% 100 == 0) {
        message("processed: ", i)
    }
}

for (i in 1:nrow(authorCitations)) {
    vidx <- which(authorCitations[i,] > 0)
    if (length(vidx) > 0 ) {
        for (v in 1:length(vidx)) {
            write(
                paste(i, "\t", vidx[v], "\t\t", authorCitations[i, vidx[v]], "\t", authorList[i,]), 
                append=T, file="data/networkcitation-2.txt")
   
        }
    }
}

stop()


# iterate over the rows (authors) and find co-authors

for (i in 1:nrow(authorPaperBiAdj)) {
    # paper indexes
    pidx <- which(authorPaperBiAdj[i,] == 1)
    # we need to sum the total citations between authors across papers
    coauthor_total_citations <- list()
    for (p in 1:length(pidx)) {
        # citations between current author and coauthors - remove the current
        # author from the list since we are interested in the edge weight
        coauthors = setdiff( which(authorPaperBiAdj[,pidx[p]] == 1), i )
        # for each coauthor count the shared citations and this becomes the 
        # edge weight between the current author and coauthor

        if (length(coauthors) > 0) {
            for (ca in 1:length(coauthors)) {
                # sum the shared citations associated with this coauthor
                citations <- sum(paperCitAdj[pidx[p],])
                akey = paste(coauthors[ca])
                if (!( akey %in% names(coauthor_total_citations))) {
                    coauthor_total_citations[[ akey ]] = 0
                }
                coauthor_total_citations[[ akey ]] = coauthor_total_citations[[ akey ]] + citations
                #message(i, "\t", akey, "\t", coauthor_total_citations[[ akey ]])
                #write( 
                #    paste(i, coauthors[ca], "", citations, sep="\t"), stdout())
            }
        }    
    }
    
    for (ca in names(coauthor_total_citations)) {
        write(
              paste(i, "\t", ca, "\t\t", coauthor_total_citations[[ ca ]]), 
              append=T, file="data/NetworkCitation.txt")
    }
}


#m <- as.matrix(paperCitAdj)
#mode(m) <- "numeric"
#g1 <- graph.adjacency(m, mode="directed")
