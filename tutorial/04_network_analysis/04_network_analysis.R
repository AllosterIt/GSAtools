######################################################################
#
#       process data
#
######################################################################
source('../../Rscripts/modules/FragMINetwork.R')

######################################################################
#       general parameters 
fragLength = 4
FDRcutoff = 0.001
percCutoff = 0.0

######################################################################
#       graphical parameters 
vertexColor = "#55555555"
vertexFrameColor = "#000000ff"
edgeColor = "#55555522"

altVertexLabelColor = "#ffffffff"
altVertexColor = "#0000ffff"
altVertexFrameColor = "#0000ffff"
altEdgeColor = "#ff00ff77"

nColor = 50

######################################################################
#
#       parameters
#
#       The script can easily be adapted to process user-provided data
#       by simply modifying the values of these parameters.
#
wSize = 1200
hSize = 800

pdbFilename = "T00.reference.pdb"
MIInputFilename = "T00.lf_MImat.out"
jointEntropyInputFilename = "T00.lf_jHmat.out"
estErrorInputFilename = "T00.lf_eeMImat.out"
pValueInputFilename = "T00.lf_pvalueMImat.out"

outputPrefix = "T04.LF"

######################################################################
#
#       load the data
 MILF = as.matrix(read.table(MIInputFilename))
 jHLF = as.matrix(read.table(jointEntropyInputFilename))
eMILF = as.matrix(read.table(estErrorInputFilename))
pMILF = as.matrix(read.table(pValueInputFilename))

pdbdm = calculateDistanceMatrix(pdbFilename, fragLength) 

nFrag = dim(MILF)[1]

######################################################################
#
#       network graph
enpdLFMatrix = filterMatrixContactLongrange(pdbdm, filterMIMatrix(renameMatrix(MILF), pMILF, eMILF, jHLF, FDRcutoff), TRUE)
write.table(enpdLFMatrix, file = paste(outputPrefix, '_enpdMI.mat', sep = ''), quote = F, row.names =F, col.names = F)

LFgs = returnGraph(pdbFilename, fragLength, 
                     MILF, pMILF, eMILF, jHLF,
                     FDRcutoff)

pcalayout = layout.norm(princomp(pdbdm)$scores[,c(1,2)], -1,1, -1,1)

outpcalayout = pcalayout
row.names(outpcalayout) = V(LFgs)$name
write.table(outpcalayout, file = paste(outputPrefix, '.pcalayout', sep = ''), quote = F, col.names = F)

write.graph(LFgs, file = paste(outputPrefix, '.graphml', sep = ''), format= "graphml")
write.graph(LFgs, file = paste(outputPrefix, '.gml', sep = ''), format= "gml")

reWeightVector = 1 - (E(LFgs)$weight - min(E(LFgs)$weight)) / (max(E(LFgs)$weight) - min(E(LFgs)$weight))
weightColorVector = rgb.palette(nColor)[as.integer(reWeightVector * nColor) + 1]
weightColorData = data.frame(get.edgelist(LFgs, names=TRUE), wc = weightColorVector)
write.table(weightColorData, file = paste(outputPrefix, '.colorData', sep = ''), quote = F, row.names = F, col.names = F)

system(paste('igraph2cytoscape.py ', outputPrefix, '.pcalayout ', outputPrefix, '.gml ', outputPrefix, '.colorData', sep =''))

######################################################################
#
#       eigenvector centrality graphs
greyPerc = 0.95
nRes = 121

#       check graphical support
if(capabilities("png")){
        png(file = paste(outputPrefix, "_graph_centrality.png", sep = ''), width = wSize, height = hSize)
        par(cex = 2)
}else{
        if(capabilities("jpeg")){
                jpeg(file = paste(outputPrefix, "_graph_centrality.jpeg", sep = ''), width = wSize, height = hSize)
                par(cex = 2)
        }else{
                bmp(file = paste(outputPrefix, "_graph_centrality.bmp", sep = ''), width = wSize, height = hSize)
                par(cex = 1.5)
        }
}

par(cex = 2)
par(las = 2)
plot(
        evcent(LFgs)$vector,
        xaxt = 'n',
        yaxt = 'n',
        xlab = 'Fragment',
        ylab = 'Eigenvector centrality',
        main = 'Centrality of the nodes in the network',
        pch = 16,
        type = 'n',
        lty = 1,
        col = rgb(0,0,1)
)
nRes = length(evcent(LFgs)$vector)
abline(h = seq(0, 1, 0.1), lty = 3, col = 'grey')
abline(v = seq(5,nFrag,5), lty = 3, col = 'grey')
points(
        evcent(LFgs)$vector,
        pch = 16,
        type = 'o',
        lty = 1,
        lwd = 4,
        col = rgb(0,0,1)
)
points(
        evcent(LFgs)$vector,
        col = rgb(1,0,0),
        pch = 20
)
axis(1, at = seq(5,(nFrag - 1),5), labels = seq(5,(nFrag - 1),5))
axis(2, seq(0.0,1.0,0.1))
dev.off()

write.table(evcent(LFgs)$vector, file = paste(outputPrefix, '.evcent.dat' , sep = ''), row.names = F, col.names = F, quote = F)

