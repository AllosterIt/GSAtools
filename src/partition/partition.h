/*
 *     $Id: partition.h 1399 2013-04-16 01:12:57Z apandini $
 *     Copyright (C) 2009-2013 Alessandro Pandini
 *
 *     This file is part of GSATools.
 *
 *     GSATools is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     GSATools is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with GSATools.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#if !defined PARTITION_H
#define PARTITION_H

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "label.h"
#include "object.h"
#include "value.h"
#include "../general/matrix.h"

/*___________________________________________________________________________*/
/* prototypes */

void partition(Collection *collection, int n, int k, float *breaks,
               int *breaks_idx, int maximize_flag, int score_idx);
void print_partition_table(Collection *collection, LabelList *labelList, int k,
                           int *breaks_idx, FILE *outputFile,
                           FILE *altOutputFile, int score_idx);
void write_breaks(float *breaks, int nBreaks, FILE *outputFile);

#endif
