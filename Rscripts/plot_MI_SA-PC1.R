######################################################################
#
#       R script to plot RMSD fitting data
#
#       input data can be obtained with
#       g_sa_encode option -rmsdlf (or -rmsdgf) and tabled with
#       the bashscripts/extract_MI_value_from_log.sh script
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
inputFilename = "extradata/example.MI.table"
outputFilenamePrefix = "figures/example.MI_SA-PC1"

######################################################################
#
#       load data
d = read.table(inputFilename, header = T)
MI = as.double(d$MI)
jH = as.double(d$jH)
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
        MI/jH,
        xlab = 'Fragment',
        ylab = expression(paste(MI^n)),
        xaxt = 'n',
        yaxt = 'n',
        type = 'l',
        pch = 16,
        col = 'white',
        ylim = c(0,0.25),
        main = 'Correlation among local motions and the first global mode (PC1)'
)
abline(h = seq(0,0.25,0.05), lty = 2, col = 'grey')
abline(v = seq(5,nFrag,5), lty = 3, col = 'grey')
axis(1, seq(5, length(MI), 5))
axis(2, seq(0, 0.25, 0.05))
lines(
        MI/jH,
        col = 'red',
        type = 'o',
        pch = 16
)
points(
        MI/jH,
        col = 'blue',
        pch = 20
)
dev.off()


