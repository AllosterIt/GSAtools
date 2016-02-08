/*
 *     $Id: mergesort.c 1398 2013-04-16 01:11:42Z apandini $
 *     Copyright (C) 2008-2013 Alessandro Pandini
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

#include "mergesort.h"

/*_____________________________________________________________________________*/
/** merge function */
void merge(unsigned char *input, int left, int right, int length, int midpoint, size_t esize, unsigned char *scratch, int (*compare) (const void *key1, const void *key2)){

    int i = 0;
    int l, r; /* counters for left and right subarrays */

    l = left; /* left subarray starts at array start */ 
    r = left + midpoint; /* right subarray starts at array midpoint */

    /* operation are performed by copy esize chunks of memory
     * indexing is dealt consistently by bytes */

    /* merge the subarrays together */ 
    /* because (l+r == length) no risk of going out of array boundaries */
    for(i = 0; i < length; i++){
        /* check if left array has elements */
        if (l < left + midpoint){
            /* check if right array has elements */
            if (r < right){
                /* compare element in position l and r */
                if (compare(&input[l * esize], &input[r * esize]) <= 0){
                    memcpy(&scratch[i * esize],  &input[l * esize],  esize);
                    l++;
                }else{
                    memcpy(&scratch[i * esize],  &input[r * esize],  esize);
                    r++;
                }
            }else{
                /* otherwise add element from left array directly */
                memcpy(&scratch[i * esize],  &input[l * esize],  esize);
                l++;
            }
        }else{
            /* otherwise add element from right array directly */
            memcpy(&scratch[i * esize],  &input[r * esize],  esize);
            r++;
        }
    }

    /* update input array with sorted subarray from scratch */
    for(i = left; i < right; i++)
        memcpy(&input[i * esize], &scratch[(i - left) * esize], esize);
}

/*_____________________________________________________________________________*/
/** mergeSort helper 
 * this is the function that actually perform the sorting where
 * left is the index of the leftmost element of the array
 * right is one past the index of the rightmost element */
void mergeSort(unsigned char *input, int left, int right, size_t esize, unsigned char *scratch, int (*compare) (const void *key1, const void *key2))
{
    int length = right - left;
    int midpoint = length/2;

    /* if the array includes only one element return to parent call */
    if(length == 1)
        return;

    /* if the array includes more than one element sorting is performed */

    /* recursively sort each subarray */
    mergeSort(input, left, left + midpoint, esize, scratch, compare);
    mergeSort(input, left + midpoint, right, esize, scratch, compare);

    /* merge the subarrays */
    merge(input, left, right, length, midpoint, esize, scratch, compare);
}

/*_____________________________________________________________________________*/
/** MergeSort */
void MergeSort(void *array, size_t size, size_t esize, int (*compare) (const void *key1, const void *key2))
{
    /* allocate memory for scratch array */
    unsigned char *scratch = safe_malloc(size * esize * sizeof(unsigned char));

    /* perform merge sorting */
    mergeSort(array, 0, size, esize, scratch, compare);

    /* free memory from scratch array */
    free(scratch);
}

