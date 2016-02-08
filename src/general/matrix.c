/*
 *     $Id: matrix.c 1398 2013-04-16 01:11:42Z apandini $
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

#include "matrix.h"

int **alloc_int_matrix(int **matrix, int nRow, int nCol)
{
	unsigned int i;

    matrix = (int **)safe_malloc(nRow * sizeof(int *));
    for (i = 0; i < nRow; ++ i)
        matrix[i] = (int *)safe_malloc(nCol * sizeof(int));

    return matrix;
}

void initialise_int_matrix(int **matrix, int nRow, int nCol, int init_value)
{
	unsigned int i,j;

    for (i = 0; i < nRow; ++ i)
	    for (j = 0; j < nCol; ++ j)
	        matrix[i][j] = init_value;

}

void print_int_matrix(int **matrix, int nRow, int nCol)
{
	int i,j;

	printf("\n");
	for (i = 0; i < nRow; i++) {
		for (j = 0; j < nCol; j++)
			printf(" %8d ", matrix[i][j]);
		printf("\n");
	}
}

void free_int_matrix(int **matrix, int nRow)
{
    unsigned int i;

    for (i = 0; i < nRow; ++ i)
        free(matrix[i]);

    free(matrix);
}

float **alloc_float_matrix(float **matrix, int nRow, int nCol)
{
	unsigned int i;

    matrix = (float **)safe_malloc(nRow * sizeof(float *));
    for (i = 0; i < nRow; ++ i)
        matrix[i] = (float *)safe_malloc(nCol * sizeof(float));

    return matrix;
}

void initialise_float_matrix(float **matrix, int nRow, int nCol, float init_value)
{
	unsigned int i,j;

    for (i = 0; i < nRow; ++ i)
	    for (j = 0; j < nCol; ++ j)
	        matrix[i][j] = init_value;

}

void print_float_matrix(float **matrix, int nRow, int nCol)
{
	int i,j;

	printf("\n");
	for (i = 0; i < nRow; i++) {
		for (j = 0; j < nCol; j++)
			printf(" %8.3f ", matrix[i][j]);
		printf("\n");
	}
}

void free_float_matrix(float **matrix, int nRow)
{
    unsigned int i;

    for (i = 0; i < nRow; ++ i)
        free(matrix[i]);

    free(matrix);
}

/*___________________________________________________________________________*/
/** 3D float matrix */
/** allocate */
float ***alloc_float_matrix3D(float ***float_matrix3D, int x, int y, int z)
{
	unsigned int i, j;

	float_matrix3D = (float ***)safe_malloc(x * sizeof(float **));
	for (i = 0; i < x; ++ i) {
		float_matrix3D[i] = (float **)safe_malloc(y * sizeof(float *));
		for (j = 0; j < y; ++ j)
			float_matrix3D[i][j] = (float *)safe_malloc(z * sizeof(float));
	}

	return float_matrix3D;
}

/** initialise */
void init_float_matrix3D(float ***float_matrix3D, int x, int y, int z, float val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				float_matrix3D[i][j][z] = val;
}


/** free */
void free_float_matrix3D(float ***float_matrix3D, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(float_matrix3D[i][j]);
        free(float_matrix3D[i]);
    }
    free(float_matrix3D);
}
