/*
 *     $Id: quatRMSD.c 1402 2013-04-16 01:19:10Z apandini $
 *     Copyright (c) 2005-2008 Douglas L. Theobald
 *     Copyright (c) 2008-2013 Alessandro Pandini
 *
 *     This file is part of GSATools. 
 *
 *     This code is a modification of QuatCharPoly.c from Douglas L. Theobald
 *     The original code can be found at:
 *     http://monkshood.colorado.edu/QCP
 *
 *     If you use this code in a publication, please reference: 
 *
 *        Douglas L. Theobald (2005)
 *        "Rapid calculation of RMSD using a quaternion-based characteristic
 *         polynomial."
 *        Acta Crystallographica A 61(4):478-480.
 *
 *        Horn, B. K. P. (1987).
 *        "Closed-form solution of absolute orientation using unit quaternions."
 *        J Opt Soc Am A 4(4):629Ð642.
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

/*____________________________________________________________________________*/
/*
    Please note that structure coordinates are stored in double **coords
    arrays. They are 3xN arrays, not Nx3 arrays as is also commonly
    used (where the x, y, z axes are interleaved). The difference is 
    something like this for storage of a structure with 8 atoms:

    Nx3: xyzxyzxyzxyzxyzxyzxyzxyz
    3xN: xxxxxxxxyyyyyyyyzzzzzzzz

    The functions can be easily modified, however, to accomodate any
    data format preference. Theobald chose this format because it is readily
    used in vectorized functions (SIMD, Altivec, MMX, SSE2, etc.).           
                                                                              */
/*____________________________________________________________________________*/
 
#include "quatRMSD.h"

double
**MatInit(const int rows, const int cols)
{
    int             i;
    double        **matrix = NULL;
    double         *matspace = NULL;

    matspace = (double *) calloc((rows * cols), sizeof(double));
    if (matspace == NULL)
    {
        perror("\n ERROR");
        puts("\n ERROR: Failure to allocate room for pointers");
        exit(EXIT_FAILURE);
    }

    /* allocate room for the pointers to the rows */
    matrix = (double **) malloc(rows * sizeof(double *));
    if (matrix == NULL)
    {
        perror("\n ERROR");
        puts("\n ERROR: Failure to allocate room for pointers");
        exit(EXIT_FAILURE);
    }

    /*  now 'point' the pointers */
    for (i = 0; i < rows; i++)
        matrix[i] = matspace + (i * cols);
    
    return(matrix);
}

void
MatDestroy(double **matrix)
{
    if (matrix[0] != NULL)
        free(matrix[0]);

    if (matrix != NULL)
        free(matrix);

    matrix = NULL;
}

void
PrintCoords(const double **coords, const int len)
{
    int             i;

    for (i = 0; i < len; ++i)
        printf("\n % 8.3f % 8.3f % 8.3f", coords[0][i], coords[1][i], coords[2][i]);
    putchar('\n');
}

/* A lot of register variables, but this sort of thing scales very well with
   new and improved processors */
static void
CalcQuarticCoeffs(const double **coords1, const double **coords2, const int len, double *coeff)
{
    double          Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double          Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
                    SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
                    SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
                    SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double          x1, x2, y1, y2, z1, z2;
    const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2],
                   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
    int             i;

    Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
    for (i = 0; i < len; ++i)
    {
        x1 = fx1[i];
        y1 = fy1[i];
        z1 = fz1[i];
        x2 = fx2[i];
        y2 = fy2[i];
        z2 = fz2[i];
   
        Sxx += (x1 * x2);
        Sxy += (x1 * y2);
        Sxz += (x1 * z2);

        Syx += (y1 * x2);
        Syy += (y1 * y2);
        Syz += (y1 * z2);

        Szx += (z1 * x2);
        Szy += (z1 * y2);
        Szz += (z1 * z2);  
    }

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    /* coeff[4] = 1.0; */
    /* coeff[3] = 0.0; */
    coeff[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    coeff[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz+Szx;
    SyzpSzy = Syz+Szy;
    SxypSyx = Sxy+Syx;
    SyzmSzy = Syz-Szy;
    SxzmSzx = Sxz-Szx;
    SxymSyx = Sxy-Syx;
    SxxpSyy = Sxx+Syy;
    SxxmSyy = Sxx-Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    coeff[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));
}

/* Evaluates the derivative of the Horn quartic for
   coefficients c and given x. */
double
eval_horn_quart_deriv(const double *c, const double x)
{
    return(2.0*(2.0*x*x + c[2])*x + c[1]);
}

/* Evaluates the Newton-Raphson correction for the Horn quartic.
   only 11 FLOPs */
static double
eval_horn_NR_corrxn(const double *c, const double x)
{
    double x2 = x*x;
    double b = (x2 + c[2])*x;
    double a = b + c[1];

    return((a*x + c[0])/(2.0*x2*x + b + a));
}

/* Evaluates the Horn quartic for coefficients c and given x. */
double
eval_horn_quart(const double *c, const double x)
{
    return(((x*x + c[2])*x + c[1])*x + c[0]);
}

/* Newton-Raphson root finding */
static double
QCProot(double *coeff, double guess, const double delta)
{
    int             i;
    double          oldg;

    for (i = 0; i < 50; ++i)
    {
        oldg = guess;
        /* guess -= (eval_horn_quart(coeff, guess) / eval_horn_quart_deriv(coeff, guess)); */
        guess -= eval_horn_NR_corrxn(coeff, guess);
    
        if (fabs(guess - oldg) < fabs(delta*guess))
            return(guess);
    }

    fprintf(stderr,
            "\n\n ERROR21: Newton-Raphson root-finding in \'QCProot()\' did not converge \n");

    exit(EXIT_FAILURE);
}

/* Calculate the inner product of some coordinates.
   This is the same as the squared radius of gyration without normalization
   for the number of atoms. */
static double
CoordsInnerProd(const double **coords, const int len)
{
    int             i;
    double          sum, tmpx, tmpy, tmpz;
    const double   *x = coords[0], *y = coords[1], *z = coords[2];

    sum = 0.0;
    for (i = 0; i < len; ++i)
    {
        tmpx = x[i];
        tmpy = y[i];
        tmpz = z[i];
        sum += (tmpx*tmpx + tmpy*tmpy + tmpz*tmpz);
    }

    return(sum);
}

void
CenterCoords(double **coords, const int len)
{
    int             i;
    double          xsum, ysum, zsum;
    double         *x = coords[0], *y = coords[1], *z = coords[2];

    xsum = ysum = zsum = 0.0;
    for (i = 0; i < len; ++i)
    {
        xsum += x[i];
        ysum += y[i];
        zsum += z[i];
    }

    xsum /= len;
    ysum /= len;
    zsum /= len;

    for (i = 0; i < len; ++i)
    {
        x[i] -= xsum;
        y[i] -= ysum;
        z[i] -= zsum;
    }
}

/* returns the sum of the squared deviations (sumdev^2) between the two structures
   rmsd = sqrt(sumdev^2 / atom_num)
*/
double
QuatCharPoly(const double **coords1, const double **coords2, const int len, double *coeff)
{
    double          innerprod;
    double          lambdamax;

    innerprod = CoordsInnerProd(coords1, len) + CoordsInnerProd(coords2, len);
    CalcQuarticCoeffs(coords1, coords2, len, coeff);
    lambdamax = QCProot(coeff, 0.5 * innerprod, 1e-6);

    return (innerprod - (2.0 * lambdamax));
}

double
QCP_rmsd(double **coords1, double **coords2, const int len, double *coeff)
{
    double          sumdev2, rmsd;

    CenterCoords(coords1, len);
    CenterCoords(coords2, len);

    sumdev2 = QuatCharPoly((const double **) coords1, (const double **) coords2, len, coeff);
    rmsd = sqrt(fabs(sumdev2 / len));

    return rmsd;
}
