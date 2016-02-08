/*
 *     $Id: transform_segment.c 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (C) 2008-2013 Alessandro Pandini and Jens Kleinjung
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

#include "transform_segment.h"

/*____________________________________________________________________________*/
/* RMSD of two coord vectors */
float coord_rmsd(Vec *v1, Vec *v2)
{
	return sqrt(pow(v1->x - v2->x, 2)
              + pow(v1->y - v2->y, 2)
			  + pow(v1->z - v2->z, 2));
}

/*____________________________________________________________________________*/
/** free memory used for structural superpositioning */
void kabsch_free(gsl_matrix *U, gsl_vector *t, gsl_matrix *X, gsl_matrix *Y){

    /* symmetry operators */
    gsl_matrix_free(U);
    gsl_vector_free(t);
    gsl_matrix_free(X);
    gsl_matrix_free(Y);
    
}

/*____________________________________________________________________________*/
/** allocate memory to symmetry operators for structural superpositioning */
void kabsch_alloc(gsl_matrix **ptr_U, gsl_vector **ptr_t, gsl_matrix **ptr_X, gsl_matrix **ptr_Y, int npoints){
    
    *ptr_U = gsl_matrix_alloc(3, 3); /* rotation matrix */
    *ptr_t = gsl_vector_alloc(3); /* translation vector */
    *ptr_X = gsl_matrix_alloc(npoints, 3); /* structure A */
    *ptr_Y = gsl_matrix_alloc(npoints, 3); /* structure B */

}

/*____________________________________________________________________________*/
/** transform atom coordinates */
__inline__ static void transform(Vec i_coords_B, Vec *i_coords_C, gsl_matrix *U, gsl_vector *t)
{
    /* create a temporary vector to hold transformed coordinates 
     * this allows source and target vectors to be the same */
    Vec atom_coords; 

    atom_coords.x =  i_coords_B.x * gsl_matrix_get(U,0,0) +
                     i_coords_B.y * gsl_matrix_get(U,0,1) +
                     i_coords_B.z * gsl_matrix_get(U,0,2) +
                     gsl_vector_get(t,0);

    atom_coords.y =  i_coords_B.x * gsl_matrix_get(U,1,0) +
                     i_coords_B.y * gsl_matrix_get(U,1,1) +
                     i_coords_B.z * gsl_matrix_get(U,1,2) +
                     gsl_vector_get(t,1);

    atom_coords.z =  i_coords_B.x * gsl_matrix_get(U,2,0) +
                     i_coords_B.y * gsl_matrix_get(U,2,1) +
                     i_coords_B.z * gsl_matrix_get(U,2,2) +
                     gsl_vector_get(t,2);

    /* update coordinates of target vector */
    i_coords_C->x =  atom_coords.x;
    i_coords_C->y =  atom_coords.y;
    i_coords_C->z =  atom_coords.z;
}

/*____________________________________________________________________________*/
/** calculate square diff of atom pairs */
__inline__ static float calc_sqdiff(Vec vector0, Vec vector1)
{
    return (  pow(vector0.x - vector1.x, 2)
            + pow(vector0.y - vector1.y, 2) 
            + pow(vector0.z - vector1.z, 2));
}

/*____________________________________________________________________________*/
/** define fragment points to superimpose */
void define_points(gsl_matrix *X, gsl_matrix *Y, Vec *coords_A, Vec *coords_B, int npoints)
{
    int i;

    for (i = 0; i < npoints; ++ i)
    {
        /** matrix Y contains points to move to = template */
        gsl_matrix_set(Y, i, 0, coords_A[i].x);
        gsl_matrix_set(Y, i, 1, coords_A[i].y);
        gsl_matrix_set(Y, i, 2, coords_A[i].z);

        /** matrix X contains points to be moved = query */
        gsl_matrix_set(X, i, 0, coords_B[i].x);
        gsl_matrix_set(X, i, 1, coords_B[i].y);
        gsl_matrix_set(X, i, 2, coords_B[i].z);
    }

}

/*____________________________________________________________________________*/
/* calculate weighted sigma square */
float sigma_square_average(gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints, int iPos, int lStruct)
{
    unsigned int i;
    float sigmaSqAve = 0.0;
    int pos;
    int wi;

    wi = npoints;

    for (i = 0; i < npoints; ++ i) 
    {
        /* weight is the number of overlapping atoms that are fitted per position */
        pos = iPos + i;
        if (pos < (npoints - 1))
            wi = pos + 1;
        else
            wi = npoints;
        if (pos  > (lStruct - npoints))
            wi = lStruct - pos;
        /* transform coordinates */
        transform(coords_B[i], &coords_C[i], U, t);
         /* calculate weighted diff2 sum */
        sigmaSqAve += (calc_sqdiff(coords_A[i], coords_C[i]) / wi);
    }

    return(sigmaSqAve);
}

/*____________________________________________________________________________*/
/** transform selected coordinates of query fragment : x' = Ux + t
    where 'U' is the rotation matrix and 't' is the translation vector
    the values of which were determined in the Kabsch routine */
float transform_coordinates(gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints)
{
    unsigned int i;
    float sqdiff_sum = 0.;

    for (i = 0; i < npoints; ++ i) 
    {
        transform(coords_B[i], &coords_C[i], U, t); /* transform coordinates */
        sqdiff_sum += calc_sqdiff(coords_A[i], coords_C[i]); /* calculate diff2 sum */
    }

    return sqrt(sqdiff_sum / (float) npoints);
}

/*____________________________________________________________________________*/
/** superimpose segment */
float superimpose_segment(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_C, int npoints)
{
    float rmsd = 0.;
    double *s = 0; /* the optimal scaling, if set != 0 */

    /** define points to match in superpositioning */
    define_points(X, Y, coords_A, coords_B, npoints);
    /** superimpose using Kabsch's algorithm */
    kabsch(npoints, X, Y, U, t, s);
    /** transform structures with the determined 'U' and 't' operators */
    rmsd = transform_coordinates(U, t, coords_A, coords_B, coords_C, npoints);

    return rmsd;
}

/*____________________________________________________________________________*/
/** extend segment */
float grow_chain(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *refcoords, int npoints, int mpoints, int overlap)
{
    int i;
    int shift_A; /* shift index for segment A */
    float rmsd = 999.;
    int extended_length;
    double *s = 0; /* the optimal scaling, if set != 0 */
    Vec *coords_C; /* temp vector to hold transformed coords */

    /** symmetry operators for structural superpositioning */
    gsl_matrix *locU = 0; /* rotation matrix */
    gsl_vector *loct = 0; /* translation vector */
    gsl_matrix *locX = 0; /* structure fragment matrix query */
    gsl_matrix *locY = 0; /* structure fragment matrix template */

    /** set final fragment lenght */
    extended_length = npoints + mpoints - overlap;

    /** allocate memory for temp vector  */
    coords_C = safe_malloc(mpoints * sizeof(Vec));
    /** this will be added to point to the C-term overlapping region of A */
    shift_A = npoints - overlap;

    /** define points to match in superpositioning */
    define_points(X, Y, coords_A + shift_A, coords_B, overlap);
    /** superimpose using Kabsch's algorithm */
    kabsch(overlap, X, Y, U, t, s);
    /** transform structures with the determined 'U' and 't' operators */
    /** rmsd is calculated against orginal position and not used further on */
    rmsd = transform_coordinates(U, t, coords_B, coords_B, coords_C, mpoints);

    /** extend C-term coords from rotated C */
    for (i = 0; i < (mpoints - overlap); i++){
        coords_A[i + npoints].x = coords_C[i + overlap].x;
        coords_A[i + npoints].y = coords_C[i + overlap].y;
        coords_A[i + npoints].z = coords_C[i + overlap].z;
    }

    /** if refer coords vector is provided, calculate rmds */
    if (refcoords != NULL){
        kabsch_alloc(&locU, &loct, &locX, &locY, extended_length);
        /** define points to match in superpositioning */
        define_points(locX, locY, refcoords, coords_A, extended_length);
        /** superimpose using Kabsch's algorithm */
        kabsch(extended_length, locX, locY, locU, loct, s);
        /** transform structures with the determined 'U' and 't' operators */
        rmsd = transform_coordinates(locU, loct, refcoords, coords_A, coords_A, extended_length);
        kabsch_free(locU, loct, locX, locY);
    }

    /** free memory from temp vector  */
    free(coords_C);

    return rmsd;
}

/*____________________________________________________________________________*/
/** extend segment */
float fitgrow_chain(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *refcoords, int npoints, int mpoints, int overlap)
{
    int i;
    int shift_A; /* shift index for segment A */
    float rmsd = 999.;
    int extended_length;
    double *s = 0; /* the optimal scaling, if set != 0 */
    Vec *temp; /* temp vector for target coords */
    Vec *coords_C; /* temp vector to hold transformed coords */

    /** symmetry operators for structural superpositioning */
    gsl_matrix *locU = 0; /* rotation matrix */
    gsl_vector *loct = 0; /* translation vector */
    gsl_matrix *locX = 0; /* structure fragment matrix query */
    gsl_matrix *locY = 0; /* structure fragment matrix template */

    /** set final fragment lenght */
    extended_length = npoints + mpoints - overlap;

    /** allocate memory for temp vectors  */
    temp = safe_malloc(mpoints * sizeof(Vec));
    coords_C = safe_malloc(mpoints * sizeof(Vec));
    /** this will be added to point to the C-term overlapping region of A */
    shift_A = npoints - overlap;

    /** fill temp vector for fitting */
    for (i = 0; i < mpoints; i++){
        if (i % 2 == 0) {
            /* even coords come from previous reconstructed trace */
            temp[i].x = coords_A[i + shift_A].x;
            temp[i].y = coords_A[i + shift_A].y;
            temp[i].z = coords_A[i + shift_A].z;
        }else{
            /* odd coords come from reference structure */
            temp[i].x = refcoords[i + shift_A].x;
            temp[i].y = refcoords[i + shift_A].y;
            temp[i].z = refcoords[i + shift_A].z;
        }
    }

    /** define points to match in superpositioning */
    define_points(X, Y, temp, coords_B, mpoints);
    /** superimpose using Kabsch's algorithm */
    kabsch(overlap, X, Y, U, t, s);
    /** transform structures with the determined 'U' and 't' operators */
    /** rmsd is calculated against orginal position and not used further on */
    rmsd = transform_coordinates(U, t, coords_B, coords_B, coords_C, mpoints);

    /** extend C-term coords from rotated C */
    for (i = 0; i < (mpoints - overlap); i++){
        coords_A[i + npoints].x = coords_C[i + overlap].x;
        coords_A[i + npoints].y = coords_C[i + overlap].y;
        coords_A[i + npoints].z = coords_C[i + overlap].z;
    }

    kabsch_alloc(&locU, &loct, &locX, &locY, extended_length);
    /** define points to match in superpositioning */
    define_points(locX, locY, refcoords, coords_A, extended_length);
    /** superimpose using Kabsch's algorithm */
    kabsch(extended_length, locX, locY, locU, loct, s);
    /** transform structures with the determined 'U' and 't' operators */
    rmsd = transform_coordinates(locU, loct, refcoords, coords_A, coords_A, extended_length);
    kabsch_free(locU, loct, locX, locY);

    /** free memory from temps vector  */
    free(temp);
    free(coords_C);

    return rmsd;
}

/*____________________________________________________________________________*/
/** extend segment */
float extend_segment(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *U, gsl_vector *t, Vec *coords_A, Vec *coords_B, Vec *coords_D, Vec *refcoords, int npoints, int mpoints, int overlap)
{
    int i;
    int shift_A; /* shift index for segment A */
    float rmsd = 999.;
    int extended_length;
    double *s = 0; /* the optimal scaling, if set != 0 */
    Vec *coords_C; /* temp vector to hold transformed coords */

    /** symmetry operators for structural superpositioning */
    gsl_matrix *locU = 0; /* rotation matrix */
    gsl_vector *loct = 0; /* translation vector */
    gsl_matrix *locX = 0; /* structure fragment matrix query */
    gsl_matrix *locY = 0; /* structure fragment matrix template */

    /** set final fragment lenght */
    extended_length = npoints + mpoints - overlap;

    /** allocate memory for temp vector  */
    coords_C = safe_malloc(mpoints * sizeof(Vec));
    /** this will be added to point to the C-term overlapping region of A */
    shift_A = npoints - overlap;

    /** define points to match in superpositioning */
    define_points(X, Y, coords_A + shift_A, coords_B, overlap);
    /** superimpose using Kabsch's algorithm */
    kabsch(overlap, X, Y, U, t, s);
    /** transform structures with the determined 'U' and 't' operators */
    /** rmsd is calculated against orginal position and not used further on */
    rmsd = transform_coordinates(U, t, coords_B, coords_B, coords_C, mpoints);

    /** copy N-term coords from A */
    for (i = 0; i < npoints; ++i ){
        coords_D[i].x = coords_A[i].x;
        coords_D[i].y = coords_A[i].y;
        coords_D[i].z = coords_A[i].z;
    }

    /** extend C-term coords from rotated C */
    for (i = 0; i < (mpoints - overlap); i++){
        coords_D[i + npoints].x = coords_C[i + overlap].x;
        coords_D[i + npoints].y = coords_C[i + overlap].y;
        coords_D[i + npoints].z = coords_C[i + overlap].z;
    }

    /** if refer coords vector is provided, calculate rmds */
    if (refcoords != NULL){
        kabsch_alloc(&locU, &loct, &locX, &locY, extended_length);
        /** define points to match in superpositioning */
        define_points(locX, locY, refcoords, coords_D, extended_length);
        /** superimpose using Kabsch's algorithm */
        kabsch(extended_length, locX, locY, locU, loct, s);
        /** transform structures with the determined 'U' and 't' operators */
        rmsd = transform_coordinates(locU, loct, refcoords, coords_D, coords_D, extended_length);
        kabsch_free(locU, loct, locX, locY);
    }

    /** free memory from temp vector  */
    free(coords_C);

    return rmsd;
}
