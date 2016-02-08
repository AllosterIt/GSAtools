/*
 *     $Id: object.c 1399 2013-04-16 01:12:57Z apandini $
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

#include "object.h"

/*____________________________________________________________________________*/
/** read the number of line in the input file */
int read_number_lines(FILE *collectionFile)
{
    int nLines = 0;
    char line[256];

    while(fgets(line, 256, collectionFile) != NULL) {
        nLines ++;
    }

    return(nLines);
}

/*____________________________________________________________________________*/
/** read the collection from the input file */
void read_collection_file(FILE *collectionFile, int n, Collection *collection)
{
    char line[256];

    collection->object[0].value = 0;
    collection->nObjects ++;

    rewind(collectionFile);

    while(fgets(line, 256, collectionFile) != NULL ) {
        collection->object[collection->nObjects].value = atof(line);
        collection->nObjects ++;
    }
}

/*____________________________________________________________________________*/
/** compare objects */
int compare_objects(Object *objectA, Object *objectB) {

    float value_A = (float) objectA->value;
    float value_B = (float) objectB->value;

    return ((value_A > value_B) ? 1 : (value_A < value_B) ? -1 : 0);
}

/*____________________________________________________________________________*/
/** calculate prefix cost */
float objectCost(Collection *collection, int i, int score_idx) {

    return(objectCost_wrapper(collection, i, score_idx));

}

/*____________________________________________________________________________*/
/** calculate partition cost */
float partitionCost(Collection *collection, int i, int x, int score_idx) {

    return(partitionCost_wrapper(collection, i, x, score_idx));

}

/*____________________________________________________________________________*/
/* initialize prefix sums */
void initialize_prefix_sums(float *p, int n, Collection *collection,
                            int score_idx) {

    initialize_prefix_sums_wrapper(p, n, collection, score_idx);

}

/*____________________________________________________________________________*/
/** calculate cost */
float costFunction(float **values, float *p, int i, int j, int x,
                   Collection *collection, int score_idx) {

    return(costFunction_wrapper(values, p, i, j, x, collection, score_idx));
}

/*____________________________________________________________________________*/
/** compare two cost values */
int compareCost(float oldValue, float newValue, int maximize_flag) {

    if (maximize_flag) {
        /* if maximization is required return 1 when newer is greater */
        if (oldValue < newValue)
            return 1;
        else
            return 0;
    } else {
        /* if minimization is required return 1 when newer is equal
                 * or smaller
         * the choice of including = in comparison differs from
         * Steve Skiena's implementation */
        if (oldValue >= newValue)
            return 1;
        else
            return 0;
    }
}

/*____________________________________________________________________________*/
/** output partition to file */
void write_partitions(Collection *collection, int start, int end,
                      FILE *outputFile)
{
    int i;

    fprintf(outputFile, "\{");
    for (i = start; i <= end; i++)
        fprintf(outputFile, " %8.3f ", collection->object[i].value);
    fprintf(outputFile, "}\n");
}
