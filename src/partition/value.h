/*
 *     $Id: value.h 1399 2013-04-16 01:12:57Z apandini $
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

#ifndef VALUES_H
#define VALUES_H

#include <stdio.h>
#include <stdlib.h>

#include "../general/safe.h"
#include "../general/mergesort.h"
#include "../statistics/probability.h"

/*____________________________________________________________________________*/
/* data structures */

typedef struct {
    float value; /* value */
    int code; /* reference code */
    int altCode; /* alternate code */
} Instance;

typedef struct {
    Instance *instance; /* input values */
    int nInstances; /* number of values */
    Set *codeSet; /* pointer to code set if available */
} Data;

/*____________________________________________________________________________*/
/* prototypes */

void read_data_file(FILE *dataFile, int n, Data *data, Set *codeSet);
void write_data(Data *data, FILE *outputFile);
void getAltCode(float *rightSideVec, int rightSideVecLength, Data *data,
                Set *altCodeSet);
void populate_probability_matrix(ProbMatrix *probMat, Data *data);
float calculate_MI_contribution(float startValue, float endValue, Data *data,
                                Set *codeSet);

#endif
