######################################################################
#
#       R script to plot transition matrix
#
#       input data can be obtained with
#       g_sa_analyze option -trans and -trmat
#
#
######################################################################
#
#       parameters
#
#       The script can easily be adapted to process user-provided data
#       by simply modifying the values of these parameters.
#
nLetters = 25
wSize = 1600
hSize = 1600
labelCutoff = 0.01
inputFilename = "../tests/data/test.lf_nMImat.out"
outputFilenamePrefix = "figures/lf_nMImat"

######################################################################
#
#       load data
d = as.matrix(read.table(inputFilename))
nFrag = dim(d)[1]
dd = d
diag(dd) = 0

#       check graphical support
if(capabilities("png")){
        png(file = paste(outputFilenamePrefix, ".png", sep = ''), width = wSize, height = hSize)
        par(cex = 2.5)
}else{
        if(capabilities("jpeg")){
                jpeg(file = paste(outputFilenamePrefix, ".jpeg", sep = ''), width = wSize, height = hSize)
                par(cex = 2.5)
        }else{
                bmp(file = paste(outputFilenamePrefix, ".bmp", sep = ''), width = wSize, height = hSize)
                par(cex = 2.0)
        }
}

#       plot data
par(mar = c(5.1, 4.1, 4.1, 6))
par(las = 2)
par(pty = 's')
par(xpd = TRUE)
image(
        dd,
        xaxt = 'n',
        yaxt = 'n',
        col = rgb(0,1,0,c(0:50)/50),
        xlab = 'Fragment',
        ylab = 'Fragment',
        main = 'normalized MI matrix'
)
axis(1, at = seq(4,(nFrag - 1),5)/(nFrag - 1), labels = seq(5,nFrag,5))
axis(2, at = seq(4,(nFrag - 1),5)/(nFrag - 1), labels = seq(5,nFrag,5))
colvals = seq(0,max(dd),max(dd)/50)
collabs = format(colvals, digits = 1)
for(i in c(0:50)){
    rect(1.05, i/50, 1.08, (i+1)/50, col = rgb(0,1,0,c(0:50)/50)[i])
    text(1.11, (i+0.5)/50, collabs[i+1], cex = 0.8)
}
dev.off()


