/*
 *     $Id: probability.h 1401 2013-04-16 01:17:37Z apandini $
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

#ifndef PROBABILITY_H_
#define PROBABILITY_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../general/matrix.h"

/*____________________________________________________________________________*/
/* data structures */

typedef struct {
    int code; /* code */
    float prob; /* probability */
} Element;

typedef struct {
    Element *element; /* element array */
    int nElements; /* number of elements */
} Set;

typedef struct {
    Set *codeSet; /* code */
    Set *altCodeSet; /* alternate code */
    float **prob; /* probability matrix */
} ProbMatrix;

/*____________________________________________________________________________*/
/* prototypes */

int compare_elements(Element *elementA, Element *elementB);
void write_probabilities(Set *codeSet, FILE *outputFile);
float Shannon(Set *codeSet);
void initialize_probability_matrix(ProbMatrix *probMat, Set *codeSet,
                                   Set *altCodeSet);
void free_probability_matrix(ProbMatrix *probMat);
float mutual_information(ProbMatrix *probMat, Set *codeSet, Set *altCodeSet);
float joint_entropy(ProbMatrix *probMat, Set *codeSet, Set *altCodeSet);

#endif /* PROBABILITY_H_ */
