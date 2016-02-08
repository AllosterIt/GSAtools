#!/bin/bash
######################################################################
#
#       check programs have been compiled
if [[ -e ../../src/g_sa_encode ]]; then
        echo "|"
else
        echo
        echo "WARNINIG: the GSATools have not been compiled."
        echo "Please read the INSTALL file for instructions."
        echo
        exit 1
fi

if [[ -e ../../src/g_sa_analyze ]]; then
        echo "|"
else
        echo
        echo "WARNINIG: the GSATools have not been compiled."
        echo "Please read the INSTALL file for instructions."
        echo 
        exit 1
fi
