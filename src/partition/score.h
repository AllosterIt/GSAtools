/*
 *     $Id: score.h 1399 2013-04-16 01:12:57Z apandini $
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

#ifndef SCORE_H
#define SCORE_H

#include <stdio.h>
#include <stdlib.h>

#include "value.h"

/*____________________________________________________________________________*/
/* data structures */

typedef struct {
    float value; /* object value */
} Object;

typedef struct {
    Object *object; /* input objects to partition */
    int nObjects; /* number of objects */
    Data *data; /* pointer to data if available */
} Collection;

/*____________________________________________________________________________*/
/* prototypes */
//float max(float a, float b);
float objectCost_wrapper(Collection *collection, int i, int score_idx);
float partitionCost_wrapper(Collection *collection, int i, int x,
                            int score_idx);
void initialize_prefix_sums_wrapper(float *p, int n, Collection *collection,
                                    int score_idx);
float costFunction_wrapper(float **values, float *p, int i, int j, int x,
                           Collection *collection, int score_idx);

#endif
