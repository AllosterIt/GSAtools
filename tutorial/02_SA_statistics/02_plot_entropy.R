######################################################################
#
#       R script to plot entropy data
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
inputFilename = "T02.lf_entropy.xvg"
outputFilenamePrefix = "T02.lf_entropy"

######################################################################
#
#       load data
d = as.matrix(read.table(inputFilename, skip = 12))
nFrag = dim(d)[1]

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
plot(
        d,
        type = 'n',
        xlab = 'Fragment',
        ylab = 'H / bit',
        xaxt = 'n',
        yaxt = 'n',
        ylim = c(0, log(25)/log(2)),
        main = 'Shannon Entropy'
)
abline(h = seq(0,log(25)/log(2),0.25), lty = 3, col = 'grey')
abline(h = seq(0,log(25)/log(2),0.5), lty = 2, col = 'grey')
abline(v = seq(5,(nFrag - 1),5), lty = 3, col = 'grey')
axis(2, seq(0,log(25)/log(2),0.5))
axis(1, at = seq(5,(nFrag - 1),5), labels = seq(5,(nFrag - 1),5))
points(
        d,
        type = 'o',
        pch = 16,
        col = 'blue'
)
points(
        d,
        pch = 20,
        col = 'red'
)
dev.off()


