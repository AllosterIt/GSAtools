######################################################################
#
#       R script to plot RMSD fitting data
#
#       
#       
#
#
######################################################################
#
#       parameters
#
#       The script can easily be adapted to process user-provided data
#       by simply modifying the values of these parameters.
#
wSize = 1200
hSize = 800
inputFilename = "T01.lf_rmsd.xvg"
outputFilenamePrefix = "T01.lf_rmsd"

######################################################################
#
#       load data
d = as.matrix(read.table(inputFilename, skip = 12), stringsAsFactors = F)
nFrag = 121
dim(d) = c(nFrag + 1,dim(d)[1]/(nFrag + 1))
# remove & character
d = d[c(1:nFrag),]
# convert to numbers
d = apply(d, c(1,2), as.double)

#       check graphical support
if(capabilities("png")){
        png(file = paste(outputFilenamePrefix, ".png", sep = ''), width = wSize, height = hSize)
        par(cex = 2)
}else{
        if(capabilities("jpeg")){
                jpeg(file = paste(outputFilenamePrefix, ".jpeg", sep = ''), width = wSize, height = hSize)
                par(cex = 2)
        }else{
                bmp(file = paste(outputFilenamePrefix, ".bmp", sep = ''), width = wSize, height = hSize)
                par(cex = 1.5)
        }
}


#       plot data
par(las = 2)
boxplot(
        t(d),
        type = 'n',
        xlab = 'Fragment',
        ylab = expression(paste("RMSD / ", ring(A))),
        xaxt = 'n',
        yaxt = 'n',
        outline = F,
        col = 'grey90',
        main = 'Local fit error (RMSD)'
)
abline(h = seq(0,max(d),0.1), lty = 2, col = 'grey')
abline(v = seq(5,(nFrag - 1),5), lty = 3, col = 'grey')
axis(2, seq(0,max(d),0.1))
axis(1, at = seq(5,(nFrag - 1),5), labels = seq(5,(nFrag - 1),5))
boxplot(
        t(d),
        outline = F,
        col = 'grey90',
        xaxt = 'n',
        yaxt = 'n',
        add = T
)
dev.off()


