/*
 *     $Id: probability.c 1401 2013-04-16 01:17:37Z apandini $
 *     Copyright (C) 2010-2013 Alessandro Pandini
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

#include "probability.h"

/*____________________________________________________________________________*/
/** compare set elements */
int compare_elements(Element *elementA, Element *elementB) {

    float value_A = (float) elementA->code;
    float value_B = (float) elementB->code;

    return ((value_A > value_B) ? 1 : (value_A < value_B) ? -1 : 0);
}

/*____________________________________________________________________________*/
/** print probability */
void write_probabilities(Set *codeSet, FILE *outputFile) {

    int i; /* index */

    for(i = 0; i < codeSet->nElements; ++i) {
        fprintf(outputFile, "%3d %8.3f\n", codeSet->element[i].code,
                codeSet->element[i].prob);
    }

}

/*____________________________________________________________________________*/
/** Shannon Entropy */
float Shannon(Set *codeSet) {

    int i; /* index */
    float H = 0.0;

    for(i = 0; i < codeSet->nElements; ++i) {
        if (codeSet->element[i].prob != 0)
            H -= codeSet->element[i].prob *
                 log(codeSet->element[i].prob) / log(2);
    }

    return(H > 0 ? H : 0);
}

/*____________________________________________________________________________*/
/** Generate Probability Matrix */
void initialize_probability_matrix(ProbMatrix *probMat, Set *codeSet,
                                   Set *altCodeSet) {

    /*________________________________________________________________________*/
    /** link code and alternate code set */
    probMat->codeSet = codeSet;
    probMat->altCodeSet = altCodeSet;

    /*________________________________________________________________________*/
    /** allocate memory for probability matrix */
    probMat->prob = alloc_float_matrix(probMat->prob,
                                       probMat->codeSet->nElements, probMat->altCodeSet->nElements);
    /*________________________________________________________________________*/
    /* initialize values and dividers matrices */
    initialise_float_matrix(probMat->prob, probMat->codeSet->nElements,
                            probMat->altCodeSet->nElements, 0.0);

}

/*____________________________________________________________________________*/
/** Free Probability Matrix */
void free_probability_matrix(ProbMatrix *probMat) {

    /*________________________________________________________________________*/
    /* free memory from matrix */
    free_float_matrix(probMat->prob, probMat->codeSet->nElements);

}

/*____________________________________________________________________________*/
/** Mutual Information */
float mutual_information(ProbMatrix *probMat, Set *codeSet, Set *altCodeSet) {

    int i,j; /* index */
    float I = 0.0;

    for(i = 0; i < codeSet->nElements; ++i)
        for(j = 0; j < altCodeSet->nElements; ++j)
            if (probMat->prob[i][j] != 0)
                I += probMat->prob[i][j] * log( (probMat->prob[i][j]) /
                                                (codeSet->element[i].prob * altCodeSet->element[j].prob) ) / log(2);

    return(I > 0 ? I : 0);
}

/*____________________________________________________________________________*/
/** Joint Entropy */
float joint_entropy(ProbMatrix *probMat, Set *codeSet, Set *altCodeSet) {

    int i,j; /* index */
    float H = 0.0;

    for(i = 0; i < codeSet->nElements; ++i)
        for(j = 0; j < altCodeSet->nElements; ++j)
            if (probMat->prob[i][j] != 0)
                H -= probMat->prob[i][j] * log(probMat->prob[i][j]) / log(2);

    return(H > 0 ? H : 0);
}
