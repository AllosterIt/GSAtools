/*
 *     $Id: label.h 1399 2013-04-16 01:12:57Z apandini $
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

#ifndef LABEL_H
#define LABEL_H

#include <stdio.h>
#include <stdlib.h>

#include "object.h"

/*____________________________________________________________________________*/
/* data structures */

typedef struct {
    char text[64]; /* label text */
} Label;

typedef struct {
    Label *label; /* label */
    int nLabels; /* number of labels */
} LabelList;

/*____________________________________________________________________________*/
/* prototypes */

void read_label_file(FILE *labelFile, int n, LabelList *labelList);

#endif
