/*
 *     $Id: matrix.h 1398 2013-04-16 01:11:42Z apandini $
 *     Copyright (C) 2007-2013 Alessandro Pandini and Jens Kleinjung
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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

#include "safe.h"

/*____________________________________________________________________________*/
/* structures */

/*___________________________________________________________________________*/
/* prototypes */

float **alloc_float_matrix(float **matrix, int nRow, int nCol);
void initialise_float_matrix(float **matrix, int nRow, int nCol, float init_value);
void print_float_matrix(float **matrix, int nRow, int nCol);
void free_float_matrix(float **matrix, int nRow);

int **alloc_int_matrix(int **matrix, int nRow, int nCol);
void initialise_int_matrix(int **matrix, int nRow, int nCol, int init_value);
void print_int_matrix(int **matrix, int nRow, int nCol);
void free_int_matrix(int **matrix, int nRow);

float ***alloc_float_matrix3D(float ***float_matrix3D, int x, int y, int z);
void init_float_matrix3D(float ***float_matrix3D, int x, int y, int z, float val);
void free_float_matrix3D(float ***float_matrix3D, int x, int y);

#endif

