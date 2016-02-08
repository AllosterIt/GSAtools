######################################################################
#
#       R script to plot data results
#
#
#
#
#
######################################################################
scriptFileList = dir('.', pattern = '.R$')

for(filename in scriptFileList){
        if(filename != 'plot_all_results.R'){
                print(paste("running script...", filename))
                source(filename)
        }
}
