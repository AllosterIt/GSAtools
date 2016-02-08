/*
 *     $Id: value.c 1399 2013-04-16 01:12:57Z apandini $
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

#include "value.h"

/*____________________________________________________________________________*/
/** read original data from the input file */
void read_data_file(FILE *dataFile, int n, Data *data, Set *codeSet)
{

    int i; /* index */
    char line[256];
    int allocated = 5;
    int added = 0;
    int (*fcmp)() = &compare_elements;

    codeSet->element = safe_malloc(allocated * sizeof(Element));
    codeSet->nElements = 0;

    rewind(dataFile);

    /* data format
     *
     * "%f %d\n"
     *
     * float space int
     *
     * */

    while(fgets(line, 256, dataFile) != NULL ) {

        assert(
            sscanf(line, "%f %d\n",
                   &(data->instance[data->nInstances].value),
                   &(data->instance[data->nInstances].code)
                  ) == 2
        );
        data->instance[data->nInstances].altCode = 0;

        /* to check input reading... uncomment the following lines */

        /*
         *  fprintf(stderr, "|%f||%d|\n",
         *  data->instance[data->nInstances].value,
         *  data->instance[data->nInstances].code);
         */

        for(i = 0; i < codeSet->nElements; ++i ) {
            if (data->instance[data->nInstances].code ==
                    codeSet->element[i].code) {
                codeSet->element[i].prob += (1.0 / n);
                added = 1;
            }
        }

        if (added == 0) {
            codeSet->element[codeSet->nElements].code =
                data->instance[data->nInstances].code;
            codeSet->element[i].prob = (1.0 / n);
            codeSet->nElements ++;

            if (codeSet->nElements == allocated) {
                allocated += 5;
                codeSet->element = safe_realloc(codeSet->element, allocated *
                                                sizeof(Element));
            }
        }

        data->nInstances++;

        added = 0;
    }

    MergeSort(codeSet->element, codeSet->nElements, sizeof(Element), fcmp);
}

/*____________________________________________________________________________*/
/** print data */
void write_data(Data *data, FILE *outputFile)
{

    int i; /* index */

    for(i = 0; i < data->nInstances; ++i) {
        fprintf(outputFile, "%8.3f %3d %3d\n", data->instance[i].value,
                data->instance[i].code, data->instance[i].altCode);
    }
}

/*____________________________________________________________________________*/
/** update alternate code */
void getAltCode(float *rightSideVec, int rightSideVecLength, Data *data,
                Set *altCodeSet)
{

    int i,j,k; /* indexes */
    int allocated = 5;
    int added = 0;
    int (*fcmp)() = &compare_elements;

    altCodeSet->element = safe_malloc(allocated * sizeof(Element));
    altCodeSet->nElements = 0;

    for(i = 0; i < data->nInstances; ++i) {
        j = 0;
        while( (data->instance[i].value > rightSideVec[j]) &
                (j < (rightSideVecLength - 1)) )
            ++j;
        data->instance[i].altCode = j;

        for(k = 0; k < altCodeSet->nElements; ++k ) {
            if (data->instance[i].altCode == altCodeSet->element[k].code) {
                altCodeSet->element[k].prob += (1.0 / data->nInstances);
                added = 1;
            }
        }

        if (added == 0) {
            altCodeSet->element[k].code = data->instance[i].altCode;
            altCodeSet->element[k].prob = (1.0 / data->nInstances);
            altCodeSet->nElements ++;

            if (altCodeSet->nElements == allocated) {
                allocated += 5;
                altCodeSet->element = safe_realloc(altCodeSet->element,
                                                   allocated * sizeof(Element));
            }
        }

        added = 0;
    }

    MergeSort(altCodeSet->element, altCodeSet->nElements, sizeof(Element),
              fcmp);
}

/*____________________________________________________________________________*/
/** check code index */
int code_index(Set *codeSet, int code) {

    int i = 0;

    while(codeSet->element[i].code != code)
        ++i;

    return(i);

}

/*____________________________________________________________________________*/
/** Populate Probability Matrix */
void populate_probability_matrix(ProbMatrix *probMat, Data *data) {

    int i, j, k;

    for(k = 0; k < data->nInstances; ++k) {
        i = code_index(probMat->codeSet, data->instance[k].code);
        j = code_index(probMat->altCodeSet, data->instance[k].altCode);
        probMat->prob[i][j] += (1.0 / data->nInstances);
    }
}

/*____________________________________________________________________________*/
/** Calculate contribution to Mutual Information */
float calculate_MI_contribution(float startValue, float endValue, Data *data,
                                Set *codeSet) {

    int i; /* index */
    float contribute = 0.0;
    float *pij;
    float pj = 0.0;

    /* intervals are defined by rightmost values */
    /*       =====|========|============|        */
    /*           [0]      [1]          [2]       */

    /* allocate memory for pij vector */
    pij = safe_malloc(sizeof(float) * codeSet->nElements);

    /* initialize pij vector */
    for(i = 0; i < codeSet->nElements; ++i)
        pij[i] = 0.0;

    /* populate pij vector */
    for(i = 0; i < data->nInstances; ++i)
        if ( (data->instance[i].value >= startValue) &
                (data->instance[i].value < endValue) ) {
            pij[code_index(codeSet, data->instance[i].code)] += (1.0 /
                    data->nInstances);
            pj += (1.0 / data->nInstances);
        }

    /* calculate contribution */
    for(i = 0; i < codeSet->nElements; ++i)
        if (pij[i] != 0)
            contribute += pij[i] * log( pij[i] / (codeSet->element[i].prob *
                                                  pj) ) / log(2);

    /** free probability matrix */
    free(pij);

    return(contribute);

}
