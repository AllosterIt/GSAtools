/*
 *     $Id: label.c 1399 2013-04-16 01:12:57Z apandini $
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

#include "label.h"

/*____________________________________________________________________________*/
/** read the labels from the input file */
void read_label_file(FILE *labelFile, int n, LabelList *labelList)
{
    char line[256];

    strcpy(labelList->label[0].text, "");

    rewind(labelFile);

    while(fgets(line, 64, labelFile) != NULL ) {
        labelList->nLabels ++;
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';
        strncpy(labelList->label[labelList->nLabels].text, line, 64);
    }
}
