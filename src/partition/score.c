/*
 *     $Id: score.c 1399 2013-04-16 01:12:57Z apandini $
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

#include "score.h"

/*____________________________________________________________________________*/
/** max of two float  */
float max(float a, float b)
{
    return( (a > b) ? a : b );
}

/*____________________________________________________________________________*/
/** MI calculate prefix cost */
float MI_objectCost(Collection *collection, int i) {

    float cost;

    cost = calculate_MI_contribution(0.0, collection->object[i].value,
                                     collection->data,
                                     collection->data->codeSet);

    return(cost);

}

/*____________________________________________________________________________*/
/** MI calculate partition cost */
float MI_partitionCost(Collection *collection, int i, int x) {

    float cost;

    cost =  calculate_MI_contribution(collection->object[i - 1].value,
                                      collection->object[x].value,
                                      collection->data,
                                      collection->data->codeSet);

    return cost;
}


/*____________________________________________________________________________*/
/* MI initialize prefix sums */
void MI_initialize_prefix_sums(float *p, int n, Collection *collection) {

    int i;

    for (i = 1; i <= n; i++)
        p[i] = MI_objectCost(collection, i);
}

/*____________________________________________________________________________*/
/** MI calculate cost */
float MI_costFunction(float **values, float *p, int i, int j, int x,
                      Collection *collection) {

    float cost;
    float Hj = 0.0;

    Hj = calculate_MI_contribution(collection->object[x].value,
                                   collection->object[i].value,
                                   collection->data,
                                   collection->data->codeSet);

    /* cost is the sum of previous contribution to MI and last one from
         * jth bin */
    cost = values[x][j-1] + Hj;

    return cost;
}

/*____________________________________________________________________________*/
/** SUM calculate prefix cost */
float SUM_objectCost(Collection *collection, int i) {

    return(collection->object[i].value);

}

/*____________________________________________________________________________*/
/** SUM calculate partition cost */
float SUM_partitionCost(Collection *collection, int i, int x) {

    int j;
    float cost = 0.0;

    for(j = i; j <= x; j ++)
        cost += SUM_objectCost(collection, j);

    return cost;
}

/*____________________________________________________________________________*/
/* SUM initialize prefix sums */
void SUM_initialize_prefix_sums(float *p, int n, Collection *collection) {

    int i;

    for (i = 1; i <= n; i++)
        p[i] = p[i-1] + SUM_objectCost(collection, i);
}

/*____________________________________________________________________________*/
/** SUM calculate cost */
float SUM_costFunction(float **values, float *p, int i, int j, int x,
                       Collection *collection) {

    float cost;

    cost = max(values[x][j-1], p[i]-p[x]);

    return cost;
}

/*____________________________________________________________________________*/
/** wrapper calculate prefix cost */
float objectCost_wrapper(Collection *collection, int i, int score_idx) {

    float returnValue;

    switch(score_idx)
    {
    case 1:
        returnValue =  SUM_objectCost(collection, i);
        break;
    case 2:
        returnValue =  MI_objectCost(collection, i);
        break;
    default:
        fprintf(stderr, "Request score is not implemented\n");
        exit(1);
    }

    return(returnValue);

}

/*____________________________________________________________________________*/
/** wrapper calculate partition cost */
float partitionCost_wrapper(Collection *collection, int i, int x,
                            int score_idx) {

    float returnValue;

    switch(score_idx)
    {
    case 1:
        returnValue =  SUM_partitionCost(collection, i, x);
        break;
    case 2:
        returnValue =  MI_partitionCost(collection, i, x);
        break;
    default:
        fprintf(stderr, "Request score is not implemented\n");
        exit(1);
    }

    return(returnValue);

}

/*____________________________________________________________________________*/
/* wrapper initialize prefix sums */
void initialize_prefix_sums_wrapper(float *p, int n, Collection *collection,
                                    int score_idx) {

    switch(score_idx)
    {
    case 1:
        SUM_initialize_prefix_sums(p, n, collection);
        break;
    case 2:
        MI_initialize_prefix_sums(p, n, collection);
        break;
    default:
        fprintf(stderr, "Request score is not implemented\n");
        exit(1);
    }
}

/*____________________________________________________________________________*/
/** wrapper calculate cost */
float costFunction_wrapper(float **values, float *p, int i, int j, int x,
                           Collection *collection, int score_idx) {

    float returnValue;

    switch(score_idx)
    {
    case 1:
        returnValue =  SUM_costFunction(values, p, i, j, x,
                                        collection);
        break;
    case 2:
        returnValue =  MI_costFunction(values, p, i, j, x,
                                       collection);
        break;
    default:
        fprintf(stderr, "Request score is not implemented\n");
        exit(1);
    }

    return(returnValue);
}
