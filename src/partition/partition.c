/*
 *     $Id: partition.c 1399 2013-04-16 01:12:57Z apandini $
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

#include "partition.h"

/*____________________________________________________________________________*/
/* print partition table */
void print_partition_table(Collection *collection, LabelList *labelList, int k,
                           int *breaks_idx, FILE *outputFile,
                           FILE *altOutputFile, int score_idx) {

    int i, j = 0;

    fprintf(outputFile, "ipart           label    value\n");
    for(i = 1; i < collection->nObjects; i++) {
        fprintf(outputFile, "%5d %15s %8.3f\n", j + 1,
                labelList->label[i].text, collection->object[i].value);
        if (breaks_idx[j] == i)
            j++;
    }

    j = 0;

    fprintf(altOutputFile, "ipart start   end     cost\n");
    fprintf(altOutputFile, "%5d %5d %5d %8.3f\n", j + 1, 1, breaks_idx[0],
            partitionCost(collection, 1, breaks_idx[0], score_idx));
    for(j = 1; j < k; j++)
        fprintf(altOutputFile, "%5d %5d %5d %8.3f\n", j + 1,
                breaks_idx[j - 1] + 1, breaks_idx[j],
                partitionCost(collection, breaks_idx[j - 1] + 1,
                              breaks_idx[j], score_idx));

}

/*____________________________________________________________________________*/
/** reconstruct the partition list */
void reconstruct_partition(Collection *collection, int **dividers, int n, int k,
                           float *breaks, int *breaks_idx)
{
    if (k==1) {
        *breaks = collection->object[n].value;
        *breaks_idx = n;
        /* AP DEBUG code */
        /* write_partitions(collection, 1, n, stderr); */
    }
    else {
        *breaks = collection->object[n].value;
        *breaks_idx = n;
        reconstruct_partition(collection, dividers, dividers[n][k], k-1,
                              (breaks + 1), (breaks_idx + 1));
        /* AP DEBUG code */
        /* write_partitions(collection, dividers[n][k]+1, n, stderr); */
    }
}

/*____________________________________________________________________________*/
/** write the breaks list */
void write_breaks(float *breaks, int nBreaks, FILE *outputFile) {

    int i; /* index */

    for(i = 0; i < (nBreaks - 1) ; ++i)
        fprintf(outputFile, "%8.3f\n", breaks[i]);

}

/*____________________________________________________________________________*/
/** reverse the breaks index list */
void reverse_breaks_idx_list(int *breaks_idx, int nBreaks,
                             int *reverseBreaks_idx) {

    int i; /* index */

    for(i = (nBreaks - 1); i >= 0 ; --i) {
        reverseBreaks_idx[nBreaks - i - 1] = breaks_idx[i];
    }

}

/*____________________________________________________________________________*/
/** reverse the breaks list */
void reverse_breaks_list(float *breaks, int nBreaks, float *reverseBreaks) {

    int i; /* index */

    for(i = (nBreaks - 1); i >= 0 ; --i) {
        reverseBreaks[nBreaks - i - 1] = breaks[i];
    }

}

/*____________________________________________________________________________*/
/** build the partition list by DP */
void partition(Collection *collection, int n, int k, float *breaks,
               int *breaks_idx, int maximize_flag, int score_idx)
{
    float *p;
    float **values = NULL;
    int **dividers = NULL;
    float cost;
    int i,j,x;
    float *tempBreaks;
    int *tempBreaks_idx;

    /*________________________________________________________________________*/
    /** allocate for breaks list */
    tempBreaks = safe_malloc((k + 1) * sizeof(float));
    tempBreaks_idx = safe_malloc((k + 1) * sizeof(int));

    /*________________________________________________________________________*/
    /** allocate memory for prefix sums, values and dividers */
    p = (float *) safe_malloc(sizeof(float) * (n + 1));
    values = alloc_float_matrix(values, (n + 1), (k + 1));
    dividers = alloc_int_matrix(dividers, (n + 1), (k + 1));

    /*________________________________________________________________________*/
    /* calculate prefix sums */
    p[0] = 0.0;
    initialize_prefix_sums(p, n, collection, score_idx);

    /*________________________________________________________________________*/
    /* initialize values and dividers matrices */
    initialise_float_matrix(values, (n + 1), (k + 1), 0.0);
    initialise_int_matrix(dividers, (n + 1), (k + 1), 0);

    /*________________________________________________________________________*/
    /* initialize boundaries */
    for (i = 1; i <= n; i++)
        values[i][1] = p[i];
    for (j = 1; j <= k; j++)
        values[1][j] = objectCost(collection, 1, score_idx);

    /*________________________________________________________________________*/
    /* calculate value matrix by DP */
    for (i = 2; i <= n; i++) {
        for (j = 2; j <= k; j++) {
            values[i][j] = maximize_flag ? FLT_MIN : FLT_MAX;
            for (x = 1; x <= (i-1); x++) {
                cost = costFunction(values, p, i, j, x,
                                    collection, score_idx);
                if (compareCost(values[i][j], cost,
                                maximize_flag) != 0) {
                    values[i][j] = cost;
                    dividers[i][j] = x;
                }
            }
        }
    }

    /*________________________________________________________________________*/
    /* AP DEBUG code */
    /* print_float_matrix(values, (n + 1), (k + 1)); */
    /* print_int_matrix(dividers, (n + 1), (k + 1)); */

    /*________________________________________________________________________*/
    /* Backtrack partition list */
    reconstruct_partition(collection, dividers, n, k, tempBreaks,
                          tempBreaks_idx);
    /* reverse breaks list */
    reverse_breaks_list(tempBreaks, k, breaks);
    reverse_breaks_idx_list(tempBreaks_idx, k, breaks_idx);

    /*________________________________________________________________________*/
    /* free memory */

    /* prefix sums */
    free(p);

    /* values and dividers matrices */
    free_float_matrix(values, (n + 1));
    free_int_matrix(dividers, (n + 1));

    /** breaks */
    free(tempBreaks);
    free(tempBreaks_idx);

}
