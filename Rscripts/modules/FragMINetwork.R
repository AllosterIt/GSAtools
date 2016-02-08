#
#   $Id: FragMINetwork.R 1390 2013-04-16 00:42:22Z apandini $
#   
#   (C) 2010-13 Alessandro Pandini
#
######################################################################
# packages
library(igraph)
library(bio3d, warn.conflicts = FALSE)
library(qvalue)

######################################################################
# Default values
fragLength = 4
FDRcutoff = 0.001

rgb.palette = colorRampPalette(c("blue", "orange", "red"), space = "rgb")

#********************************************************************#
# Functions to:                                                      #
#               process matrices                                     #
#********************************************************************#

######################################################################
# rename row and column with R_number
renameMatrix = function(m){
        
        colnames(m) = sapply(
                c(1:dim(m)[2]), 
                function(x){return(paste('R_', x, sep = ''))}
                )
        rownames(m) = colnames(m)
        
        return(m)
}

#********************************************************************#
# Functions to:                                                      #
#               calculate distance based statistics of MI values     #
#********************************************************************#

######################################################################
# calculate distance matrix from pdb file from 1 to (n - l)
calculateDistanceMatrix = function(pdbfilename, fragLength){

        pdb = read.pdb(pdbfilename)
        pdbDm = dm(pdb, selection = "calpha")
        nFrag = dim(pdbDm)[1] - fragLength + 1
        pdbDmFrag = pdbDm[c(1:nFrag), c(1:nFrag)]

        return(pdbDmFrag)
}

#********************************************************************#
# Functions to:                                                      #
#               filter the MI matrix                                 #
#********************************************************************#

######################################################################
# return kth p-value according to Benjamini–Hochberg 
pvalueCutoff = function(pvalue, alpha){
        
        m = length(pvalue)
        pk = sort(pvalue)
        kma  = alpha * c(1:m) / m  
        
        k = sum(pk <= kma)
        
        return(pk[k])
}

######################################################################
# calculate p-value cutoff according to FDRthreshold 
filterValue = function(pvalueVec, FDRcutoff){
        
        pcutoff =  pvalueCutoff(pvalueVec, FDRcutoff)
        
        return(pcutoff)
}

######################################################################
# filter the MI matrix by significance after FDR test 
filterMIMatrix = function(MIMatrix, pMIMatrix, eMIMatrix, HMatrix, FDRcutoff){
        
        neMIMatrix = (MIMatrix - eMIMatrix) / HMatrix
        
        pVec = as.vector(pMIMatrix[upper.tri(pMIMatrix)])
        
        pcutoff = filterValue(pVec, FDRcutoff) 
        
        filteredMIMatrix = neMIMatrix
        filteredMIMatrix[pMIMatrix > pcutoff] = 0
        
        return(filteredMIMatrix)
}

######################################################################
# filter the MI matrix by contact/longrange cutoffs 
filterMatrixContactLongrange = function(pdbDmFrag, enMIMat, printFlag, LRPerc = 0.95, COPerc = 0.75){ 

        MIvec = as.vector(enMIMat[upper.tri(enMIMat)])
        Dvec  = as.vector(pdbDmFrag[upper.tri(pdbDmFrag)])
        
        MIcontact = MIvec[Dvec > 4 & Dvec <= 12]
        MIlongrange = MIvec[Dvec > 12]

        coCutoff = quantile(MIcontact, COPerc)
        lrCutoff = quantile(MIlongrange, LRPerc)

        newMat = enMIMat
        newMat[pdbDmFrag <= 4] = 0 
        newMat[pdbDmFrag > 4 & pdbDmFrag <= 12 & enMIMat <= coCutoff] = 0 
        newMat[pdbDmFrag > 12 & enMIMat <= lrCutoff] = 0 

        if(printFlag){
                print(paste('Contact cutoff    [at ', COPerc, ']: ', format(coCutoff, digits = 3), sep = ''))
                print(paste('Long range cutoff [at ', LRPerc, ']: ', format(lrCutoff, digits = 3), sep = ''))
                print(paste('non-zero contact region values:      ', sum(MIcontact > coCutoff),sep = ''))
                print(paste('non-zero long range region values:   ', sum(MIlongrange > lrCutoff),sep = ''))
                print(paste('non-zero values:                     ', sum(newMat[upper.tri(newMat)] > 0),sep = ''))
                        
        }

        return(newMat)

}

#********************************************************************#
# Functions to:                                                      #
#               build the graph representation                       #
#********************************************************************#

######################################################################
# build the graph from adjecent matrix with w = 1-I 
buildGraphWeightComplement = function(adjMatrix){
        
        g = graph.adjacency(
                        adjMatrix,
                        weighted  = T,
                        mode = 'undirected',
                        diag = F)
        E(g)$weight = 1 - E(g)$weight
        return(g)
}

######################################################################
# build the graph representation 
returnGraph = function(pdbfilename, fragLength, MIMatrix, pMIMatrix, eMIMatrix, jHMatrix, FDRcutoff, LRPerc = 0.95, COPerc = 0.75){
        
        pdbDmFrag= calculateDistanceMatrix(pdbfilename, fragLength)
        adjMIMatrix = filterMIMatrix(renameMatrix(MIMatrix), pMIMatrix, eMIMatrix, jHMatrix, FDRcutoff)
        distAdjMIMatrix = filterMatrixContactLongrange(pdbDmFrag, adjMIMatrix, FALSE, LRPerc, COPerc)
        gs = buildGraphWeightComplement(distAdjMIMatrix)
        dgs = addCAdistEdge(gs, pdbDmFrag)
        
        return(dgs)
}

#********************************************************************#
# Functions to:                                                      #
#               process edge and vertices of a graph                 #
#********************************************************************#

######################################################################
# add CA dist attributes to edges in a graph
addCAdistEdge = function(g, pdbDmFrag){

        distVector = apply(
                        get.edgelist(g),
                        1,
                        function(x){
                                      fi = strsplit(x[1], '_')[[1]][2]
                                      fj = strsplit(x[2], '_')[[1]][2]
                                      return(pdbDmFrag[fi,fj])
                              }       
                      ) 

        E(g)$cad = distVector

        return(g)
}

######################################################################
