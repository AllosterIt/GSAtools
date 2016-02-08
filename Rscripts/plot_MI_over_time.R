######################################################################
#
#       R script to plot RMSD fitting data
#
#       input data can be obtained with
#       g_sa_encode option -rmsdlf (or -rmsdgf)
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

SAColor = c(
        # A         B         C         D         E
        "#0000FF","#0035FF","#006AFF","#009FFF","#00D4FF",
        # F         G         H         I         J
        "#05FFF4","#1FFFBF","#3AFF8A","#54FF55","#6FFF1F",
        # K         L         M         N         O
        "#89FF00","#A4FF00","#BEFF00","#D9FF00","#F4FF00",
        # P         Q         R         S         T
        "#FFF300","#FFE100","#FFCE00","#FFBB00","#FFA800",
        # U         V         W         X         Y
        "#FF8900","#FF6700","#FF4400","#FF2200","#FF0000"
)

inputFilename = "../tests/data/test.lf_MI.80ns.xvg"
outputFilenamePrefix = "figures/lf_MI.PC1-Fragment80"

######################################################################
#
#       load data
d = read.table(inputFilename, skip = 21)
names(d) = c("time", "value", "SA", "partition")
nFrames = dim(d)[1]

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
par(xpd = T)
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(
        d$partition,
        type = 'n',
        xlab = 'Time / ns',
        ylab = 'Partition n.',
        xaxt = 'n',
        yaxt = 'n',
        main = 'Partition of functional value (according to the encoding of Fragment 80)'
)
axis(2, at = seq(0,max(d$partition),1), labels = seq(1,(max(d$partition) + 1),1))
axis(1, at = seq(5000,nFrames,5000), labels = seq(5,(nFrames/1000),5))
points(
        d$partition,
        pch = 16,
        col = SAColor[d$SA],
        xaxt = 'n',
        yaxt = 'n'
)
        for(i in c(0:24)){
                rect(   nFrames + 5000,
                        i,
                        nFrames + 8000,
                        (i + 1),
                        col = SAColor[i + 1], border = 'white')
                text(   nFrames + 6500,
                        (i + 0.5),
                        LETTERS[i + 1],
                        cex = 1)
        }
dev.off()


