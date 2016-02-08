#!/bin/bash
######################################################################
#
#       check programs have been compiled
if [[ -e ../src/g_sa_encode ]]; then
        echo
else
        echo
        echo "WARNINIG: the GSATools have not been compiled."
        echo "Please read the INSTALL file for instructions."
        echo
        exit
fi

if [[ -e ../src/g_sa_analyze ]]; then
        echo
else
        echo
        echo "WARNINIG: the GSATools have not been compiled."
        echo "Please read the INSTALL file for instructions."
        echo 
        exit
fi

######################################################################
#
#       set GMX environment variables
source ../scripts/GSATOOLSRC.bash

echo "--------------------------------------------------"
echo " To run the tutorial a GSATOOLSRC script has been "
echo " sourced to set the correct GROMACS environmental "
echo " variables.                                       "
echo
echo " Please remember to source the GSATOOLSRC before  " 
echo " using the g_sa_encode and/or g_sa_analyze.       " 
echo " You can find the script suitable for your shell  " 
echo " in the scripts directory.                        " 
echo

######################################################################
#
#       Tutorial
#
######################################################################
echo "---------------------------------------------------------------"
echo "|"
echo "|         01 SA encoding"
echo "|         results available in the subdir 01_SA_encoding"
echo "|"
echo "|         running:"
cd 01_SA_encoding
/bin/bash 01_SA_encoding.sh
cd ..
echo
echo "---------------------------------------------------------------"
echo "|"
echo "|         02 SA statistics"
echo "|         results available in the subdir 02_SA_statistics"
echo "|"
echo "|         running:"
cd 02_SA_statistics
/bin/bash 02_SA_statistics.sh
cd ..
echo
echo "---------------------------------------------------------------"
echo "|"
echo "|         03 local correlation"
echo "|         results available in the subdir 03_local_correlation"
echo "|"
echo "|         running:"
cd 03_local_correlation
/bin/bash 03_local_correlation.sh
cd ..
echo 
echo "---------------------------------------------------------------"
echo "|"
echo "|         04 network analysis"
echo "|         results available in the subdir 04_network_analysis"
echo "|"
echo "|         running:"
cd 04_network_analysis
/bin/bash 04_network_analysis.sh
cd ..
echo 
echo "---------------------------------------------------------------"
echo "|"
echo "|         05 functional analysis"
echo "|         results available in the subdir 05_functional_analysis"
echo "|"
echo "|         running:"
cd 05_functional_analysis
/bin/bash 05_functional_analysis.sh
cd ..
echo
echo "---------------------------------------------------------------"
echo 
echo " To run the tutorial a GSATOOLSRC script has been "
echo " sourced to set the correct GROMACS environmental "
echo " variables.                                       "
echo
echo " Please remember to source the GSATOOLSRC before  " 
echo " using the g_sa_encode and/or g_sa_analyze.       " 
echo " You can find the script suitable for your shell  " 
echo " in the scripts directory.                        " 
echo

