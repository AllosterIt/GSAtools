/*
 *     $Id: object.h 1399 2013-04-16 01:12:57Z apandini $
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

#if !defined OBJECT_H
#define OBJECT_H

#include <stdio.h>
#include <stdlib.h>

#include "score.h"
#include "value.h"

/*____________________________________________________________________________*/
/* data structures */

/*____________________________________________________________________________*/
/* prototypes */

int read_number_lines(FILE *collectionFile);
void read_collection_file(FILE *collectionFile, int n, Collection *collection);
int compare_objects(Object *objectA, Object *objectB);
float objectCost(Collection *collection, int i, int score_idx);
float partitionCost(Collection *collection, int i, int x, int score_idx);
void initialize_prefix_sums(float *p, int n, Collection *collection,
                            int score_idx);
float costFunction(float **values, float *p, int i, int j, int x,
                   Collection *collection, int score_idx);
int compareCost(float oldValue, float newValue, int maximize_flag);
void write_partitions(Collection *collection, int start, int end,
                      FILE *outputFile);

#endif
