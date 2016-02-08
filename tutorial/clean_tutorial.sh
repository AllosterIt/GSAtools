#!/bin/bash
#       remove file from stage 01
rm -f 01_SA_encoding/T01.*
rm -f 01_SA_encoding/*.Rout
#       remove file from stage 02
rm -f 02_SA_statistics/T02.*
rm -f 02_SA_statistics/*.Rout
#       remove file from stage 03
rm -f 03_local_correlation/T03.*
rm -f 03_local_correlation/*.Rout
#       remove file from stage 04
rm -f 04_network_analysis/T04.*
rm -f 04_network_analysis/*.Rout
#       remove file from stage 05
rm -f 05_functional_analysis/T05.*
rm -f 05_functional_analysis/*.Rout
