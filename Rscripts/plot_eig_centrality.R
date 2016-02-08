######################################################################
#
#       R script to plot transition matrix
#
#       input data can be obtained with
#       g_sa_analyze option -MImat -eeMImat -jHmat -pvalueMImat
#
#
######################################################################
source('modules/FragMINetwork.R')

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

pdbFilename = "extradata/example.reference.pdb"
MIInputFilename = "extradata/example.lf_MImat.out"
jointEntropyInputFilename = "extradata/example.lf_jHmat.out"
estErrorInputFilename = "extradata/example.lf_eeMImat.out"
pValueInputFilename = "extradata/example.lf_pvalueMImat.out"

outputFilenamePrefix = "figures/example.LF"

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

LFgs = returnGraph(pdbFilename, fragLength, 
                     MILF, pMILF, eMILF, jHLF,
                     FDRcutoff)

######################################################################
#
#       eigenvector centrality graph
greyPerc = 0.95
nRes = 121

#       check graphical support
if(capabilities("png")){
        png(file = paste(outputFilenamePrefix, "_graph_centrality.png", sep = ''), width = wSize, height = hSize)
        par(cex = 2)
}else{
        if(capabilities("jpeg")){
                jpeg(file = paste(outputFilenamePrefix, "_graph_centrality.jpeg", sep = ''), width = wSize, height = hSize)
                par(cex = 2)
        }else{
                bmp(file = paste(outputFilenamePrefix, "_graph_centrality.bmp", sep = ''), width = wSize, height = hSize)
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

