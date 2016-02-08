######################################################################
#
#       R script to plot profile data
#
#       input data can be obtained with
#       g_sa_encode option -proflf (or -profgf)
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
hSize = 1000
inputFilename = "../tests/data/test.lf_prof.dat"
outputFilenamePrefix = "figures/lf_prof"

######################################################################
#
#       load data
d = as.matrix(read.table(inputFilename))
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
par(mar = c(5.1, 4.1, 4.1, 6))
par(xpd = TRUE)
image(
        d,
        xaxt = 'n',
        yaxt = 'n',
        col = gray(c(50:1)/50),
        xlab = 'Fragment',
        ylab = 'SA Letter',
        main = 'Structural sequence profile'
)
dy0 = 1 / (nLetters * 2)
dx0 = 1 / (nFrag * 2)
segments(seq(0,1,1/(nFrag - 1)) + (1/(nFrag - 1) * 0.5), y0 = -dy0, y1 = 1 + dy0, lty = 3, col = 'green')
segments(x0 = -dx0, x1 = 1 + dx0, y0 = seq(0,1,1/(nLetters - 1)) + (1/(nLetters - 1) * 0.5), lty = 3, col = 'green')
axis(2, at = c(0:(nLetters - 1)) / (nLetters - 1), labels = LETTERS[c(1:nLetters)])
axis(1, at = seq(4,(nFrag - 2),5) / (nFrag - 1), labels = seq(5,(nFrag - 1),5))
collabs = format(seq(0,1,0.02), digits = 1)
for(i in c(0:50)){
    rect(1.05, i/50, 1.08, (i+1)/50, col = gray(c(50:1)/50)[i])
    text(1.11, (i+0.5)/50, collabs[i+1], cex = 0.7)
    
}
dev.off()


