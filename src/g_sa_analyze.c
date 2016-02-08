/*
 *     $Id: g_sa_analyze.c 1408 2013-04-16 23:33:30Z apandini $
 *     Copyright (C) 2011-2013 Alessandro Pandini
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

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "pdbio.h"
#include "xvgr.h"
#include "tpxio.h"
#include "matio.h"
#include "float.h"

#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_float.h>

#include "general/safe.h"
#include "partition/partition.h"
#include "sequence/sequence.h"
#include "statistics/code.h"
#include "statistics/probability.h"
#include "structure/fragments.h"
#include "structure/fragment_colour.h"
#include "structure/pdb_structure.h"
#include "structure/transform_segment.h"

/*____________________________________________________________________________*/
/* global variables */
#ifdef MPI
#include <mpi.h>
int my_rank; /* rank of 'this' node */
int nodes; /* number of nodes */
#else
int my_rank = 0;
int nodes = 1;
#endif

/*____________________________________________________________________________*/
/** calculate estimated Mutual Information error */
float estimate_MI_error(ProbMatrix *probMat, Set *codeSet, Set *altCodeSet, int nObs) {

    int i,j;
    int Mxy =0;
    int Mx = 0;
    int My = 0;
    float eeMI = 0.0;

    for(i = 0; i < probMat->codeSet->nElements; ++i)
        if (probMat->codeSet->element[i].prob > 0.0)
            Mx ++;

    for(i = 0; i < probMat->altCodeSet->nElements; ++i)
        if (probMat->altCodeSet->element[i].prob > 0.0)
            My ++;

    for(i = 0; i < probMat->codeSet->nElements; ++i)
        for(j = 0; j < probMat->altCodeSet->nElements; ++j)
            if (probMat->prob[i][j] > 0.0)
                Mxy ++;

    eeMI = ((float) (Mxy - Mx - My + 1)) / ((float) 2 * nObs );

    return(eeMI);

}

/*____________________________________________________________________________*/
/** calculate Joint entropy for two columns in fasta alignment */
float column_joint_entropy(SeqSet *fastaSequenceSet, ProbMatrix *probMat, int icol, int jcol) {

    int i,k,l; /* indexes */
    char *iString, *jString; /* string placeholders */
    float Jentropy = 0.0; /* Joint entropy value */

    /* create string from column */
    iString = get_string_from_column(fastaSequenceSet, icol);
    jString = get_string_from_column(fastaSequenceSet, jcol);

    /* initialize probability matrix to 0.0 */
    initialise_float_matrix(probMat->prob, probMat->codeSet->nElements, probMat->altCodeSet->nElements, 0.0);
    /* initialize probability vectors to 0.0 */
    for(k = 0; k < probMat->codeSet->nElements; ++k)
        probMat->codeSet->element[k].prob = 0.0;
    for(l = 0; l < probMat->altCodeSet->nElements; ++l)
        probMat->altCodeSet->element[l].prob = 0.0;

    for(i = 0; i < fastaSequenceSet->nSequences; ++i) {
        for(k = 0; k < probMat->codeSet->nElements; ++k) {
            if (probMat->codeSet->element[k].code == iString[i]) {
                probMat->codeSet->element[k].prob += (1.0 / fastaSequenceSet->nSequences);
                break;
            }
        }
        for(l = 0; l < probMat->altCodeSet->nElements; ++l) {
            if (probMat->altCodeSet->element[l].code == jString[i]) {
                probMat->altCodeSet->element[l].prob += (1.0 / fastaSequenceSet->nSequences);
                break;
            }
        }
        probMat->prob[k][l] += (1.0 / fastaSequenceSet->nSequences);
    }

    Jentropy = joint_entropy(probMat, probMat->codeSet, probMat->altCodeSet);

    free(iString);
    free(jString);

    return(Jentropy);

}

/*____________________________________________________________________________*/
/** calculate randomized Mutual Information */
void random_column_mutual_information(float *ptr_meanMI, float *ptr_stdMI, float *ptr_pValueMI, float MIij, SeqSet *fastaSequenceSet, ProbMatrix *probMat, int icol, int jcol, int nSurrogates,	gsl_rng *rndGenerator) {

    int i,k,l,m; /* indexes */
    char *iString, *jString; /* string placeholders */
    float MI = 0.0; /* Mutual Information value */
    float *MIvector; /* vector of Mutual Information values */
    float nZgreater = 0; /* number of Z score values greater than Zij */
#ifdef MPI
    int dimx = (int)ceil(nSurrogates / nodes); /* dimension of MPI loop */
    float nodeMIvector[dimx]; /* MI vector holding results of this node */
    /* allocate memory for Mutual Information vector */
    MIvector = safe_malloc((dimx * nodes) * sizeof(float));
#else
    /* allocate memory for Mutual Information vector */
    MIvector = safe_malloc(nSurrogates * sizeof(float));
#endif
    /* create string from column */
    iString = get_string_from_column(fastaSequenceSet, icol);
    jString = get_string_from_column(fastaSequenceSet, jcol);

#ifdef MPI
    for(m = 0; m < dimx; ++m) {
        /* skip excess loops */
        if (((my_rank * dimx) + m) < nSurrogates) {
#else
    for(m = 0; m < nSurrogates; ++m) {
#endif
            gsl_ran_shuffle(rndGenerator, iString, fastaSequenceSet->nSequences, sizeof (char));

            /* initialize probability matrix to 0.0 */
            initialise_float_matrix(probMat->prob, probMat->codeSet->nElements, probMat->altCodeSet->nElements, 0.0);
            /* initialize probability vectors to 0.0 */
            for(k = 0; k < probMat->codeSet->nElements; ++k)
                probMat->codeSet->element[k].prob = 0.0;
            for(l = 0; l < probMat->altCodeSet->nElements; ++l)
                probMat->altCodeSet->element[l].prob = 0.0;

            for(i = 0; i < fastaSequenceSet->nSequences; ++i) {
                for(k = 0; k < probMat->codeSet->nElements; ++k) {
                    if (probMat->codeSet->element[k].code == iString[i]) {
                        probMat->codeSet->element[k].prob += (1.0 / fastaSequenceSet->nSequences);
                        break;
                    }
                }
                for(l = 0; l < probMat->altCodeSet->nElements; ++l) {
                    if (probMat->altCodeSet->element[l].code == jString[i]) {
                        probMat->altCodeSet->element[l].prob += (1.0 / fastaSequenceSet->nSequences);
                        break;
                    }
                }
                probMat->prob[k][l] += (1.0 / fastaSequenceSet->nSequences);
            }

            MI = mutual_information(probMat, probMat->codeSet, probMat->altCodeSet);

#ifdef MPI
            nodeMIvector[m] = MI;

        }
        else {
            nodeMIvector[m] = 0.;
        }
    }

    /* communicate data between all nodes */
    MPI_Allgather(nodeMIvector, dimx, MPI_FLOAT, MIvector, dimx, MPI_FLOAT, MPI_COMM_WORLD);

    for (m = 0; m < nSurrogates; ++ m)
        if (MIvector[m] > MIij)
            nZgreater ++;
#else
            MIvector[m] = MI;

            if (MI > MIij)
                nZgreater ++;
        }
#endif

    *ptr_meanMI = gsl_stats_float_mean(MIvector, 1, nSurrogates);
    *ptr_stdMI = gsl_stats_float_sd(MIvector, 1, nSurrogates);

    *ptr_pValueMI = ((float) nZgreater / (float) nSurrogates);

    free(iString);
    free(jString);

    free(MIvector);

}

/*____________________________________________________________________________*/
/** calculate Mutual Information for two columns in fasta alignment */
float column_mutual_information(SeqSet *fastaSequenceSet, ProbMatrix *probMat, int icol, int jcol, float *ptr_eeMI) {

    int i,k,l; /* indexes */
    char *iString, *jString; /* string placeholders */
    float MI = 0.0; /* Mutual Information value */

    /* create string from column */
    iString = get_string_from_column(fastaSequenceSet, icol);
    jString = get_string_from_column(fastaSequenceSet, jcol);

    /* initialize probability matrix to 0.0 */
    initialise_float_matrix(probMat->prob, probMat->codeSet->nElements, probMat->altCodeSet->nElements, 0.0);
    /* initialize probability vectors to 0.0 */
    for(k = 0; k < probMat->codeSet->nElements; ++k)
        probMat->codeSet->element[k].prob = 0.0;
    for(l = 0; l < probMat->altCodeSet->nElements; ++l)
        probMat->altCodeSet->element[l].prob = 0.0;

    for(i = 0; i < fastaSequenceSet->nSequences; ++i) {
        for(k = 0; k < probMat->codeSet->nElements; ++k) {
            if (probMat->codeSet->element[k].code == iString[i]) {
                probMat->codeSet->element[k].prob += (1.0 / fastaSequenceSet->nSequences);
                break;
            }
        }
        for(l = 0; l < probMat->altCodeSet->nElements; ++l) {
            if (probMat->altCodeSet->element[l].code == jString[i]) {
                probMat->altCodeSet->element[l].prob += (1.0 / fastaSequenceSet->nSequences);
                break;
            }
        }
        assert((k < probMat->codeSet->nElements) && "k beyond matrix limit!");
        assert((l < probMat->altCodeSet->nElements) && "l beyond matrix limit!")
        probMat->prob[k][l] += (1.0 / fastaSequenceSet->nSequences);
    }

    MI = mutual_information(probMat, probMat->codeSet, probMat->altCodeSet);
    *ptr_eeMI = estimate_MI_error(probMat, probMat->codeSet, probMat->altCodeSet, fastaSequenceSet->nSequences);

    free(iString);
    free(jString);

    return(MI);

}

/*____________________________________________________________________________*/
/** calculate transition frequencies from input sequence set */
void calculate_transition_frequencies(SeqSet *inputSequenceSet, ProbMatrix *transition_matrix, int isize, atom_id *index) {

    int i,j,k,l,s; /* indexes */
    char *iString; /* string placeholder */
    int nTransition; /* number of transitions */
    int sequenceLength = 0; /* sequence length */
    int *index_vector = 0; /* vector of position indeces */

    /* sequence length */
    if (isize == 0) {
        sequenceLength = inputSequenceSet->sequence[0].length;
        index_vector = safe_malloc(sizeof(int) * sequenceLength);
        for(i = 0; i < sequenceLength; ++i)
            index_vector[i] = i;
    } else {
        sequenceLength = isize;
        index_vector = safe_malloc(sizeof(int) * sequenceLength);
        for(i = 0; i < sequenceLength; ++i)
            index_vector[i] = index[i];
    }

    /* number of transitions */
    nTransition =  inputSequenceSet->nSequences - 1;

    /* initialize probability vectors to 0.0 */
    for(k = 0; k < transition_matrix->codeSet->nElements; ++k)
        transition_matrix->codeSet->element[k].prob = 0.0;
    for(l = 0; l < transition_matrix->altCodeSet->nElements; ++l)
        transition_matrix->altCodeSet->element[l].prob = 0.0;

    /* calculate code set probabilities */
    for(s = 0; s < sequenceLength; ++s) {
        i = index_vector[s];
        /* create string from column */
        iString = get_string_from_column(inputSequenceSet, i);
        for(j = 0; j < nTransition; ++j) {
            for(k = 0; k < transition_matrix->codeSet->nElements; ++k) {
                if (transition_matrix->codeSet->element[k].code == iString[j])
                    transition_matrix->codeSet->element[k].prob += (1.0 / (sequenceLength * nTransition));
            }
        }
        free(iString);
    }

    for(s = 0; s < sequenceLength; ++s) {
        i = index_vector[s];
        /* create string from column */
        iString = get_string_from_column(inputSequenceSet, i);
        for(j = 0; j < nTransition; ++j) {
            for(k = 0; k < transition_matrix->codeSet->nElements; ++k)
                if (transition_matrix->codeSet->element[k].code == iString[j])
                    break;
            for(l = 0; l < transition_matrix->altCodeSet->nElements; ++l)
                if (transition_matrix->altCodeSet->element[l].code == iString[j + 1])
                    break;
            transition_matrix->prob[k][l] += (1.0 / (sequenceLength * nTransition * transition_matrix->codeSet->element[k].prob));
        }
        free(iString);
    }

    free(index_vector);
}

/*____________________________________________________________________________*/
/** initialize xpm mapping */
void initialize_xpm_SA_mapping(t_mapping *SA_mapping, FragmentSet *fragment_set, SA_colour *colours) {

    int i;
    char SAstring[2];

    for(i = 0; i < fragment_set->nFragment; ++i) {
        SAstring[0] = fragment_set->codeOrder[i];
        SAstring[1] = '\0';
        SA_mapping[i].code.c1 = fragment_set->codeOrder[i];
        SA_mapping[i].code.c2 = 0;
        SA_mapping[i].desc = strdup(SAstring);
        if (colours != NULL) {
            SA_mapping[i].rgb.r = colours->colour[i].r;
            SA_mapping[i].rgb.g = colours->colour[i].g;
            SA_mapping[i].rgb.b = colours->colour[i].b;
        } else {
            SA_mapping[i].rgb.r = ((float) (fragment_set->nFragment - (i + 1))) / fragment_set->nFragment;
            SA_mapping[i].rgb.g = 0.0;
            SA_mapping[i].rgb.b = 1.0;
        }
    }

}

/*____________________________________________________________________________*/
/** initialize xpm matrix */
void initialize_xpm_matrix(t_matrix *xpm_mat, int nx, int ny) {

    int i;

    xpm_mat->nx = nx;
    xpm_mat->ny = ny;
    xpm_mat->bDiscrete = TRUE;
    xpm_mat->axis_x = 0;
    xpm_mat->axis_y = 0;
    srenew(xpm_mat->axis_x, nx);
    for(i = 0; i < nx; ++i)
        xpm_mat->axis_x[i] = (float) i;
    srenew(xpm_mat->axis_y, ny);
    for(i = 0; i < ny; ++i)
        xpm_mat->axis_y[i] = (float) i;

}

/*____________________________________________________________________________*/
/** free xpm matrix */
void free_xpm_matrix(t_matrix *xpm_mat, int nmatcol, int nFragment) {

    int i;

    sfree(xpm_mat->axis_x);
    sfree(xpm_mat->axis_y);

    sfree(xpm_mat->map);

    for(i = 0; i < nmatcol; ++i)
        sfree(xpm_mat->matrix[i]);
    sfree(xpm_mat->matrix);

}

/*____________________________________________________________________________*/
/** create and save xpm matrix */
void save_SA_xpm_matrix(FILE *xpmout, SeqSet *fastaSequenceSet, int sequenceLength, FragmentSet *fragment_set, SA_colour *SA_colour_set, int nColourSets) {

    int i,j; /* index */
    t_matrix xpm_mat; /* matrix for xpm printout */
    t_xpmelmt c; /* placeholder for mapping instance */
    int colourIndex = -1;

    initialize_xpm_matrix(&xpm_mat, fastaSequenceSet->nSequences, sequenceLength);
    xpm_mat.map = safe_malloc(sizeof(t_mapping) * fragment_set->nFragment);

    for(i = 0; i < nColourSets; ++i)
        if (strcmp(SA_colour_set[i].alphabetName, fragment_set->setname) == 0)
            colourIndex = i;
    if (colourIndex > -1)
        initialize_xpm_SA_mapping(xpm_mat.map, fragment_set, &(SA_colour_set[colourIndex]));
    else
        initialize_xpm_SA_mapping(xpm_mat.map, fragment_set, NULL);

    xpm_mat.nmap = fragment_set->nFragment;

    sprintf(xpm_mat.title, "SA xpm Matrix");
    sprintf(xpm_mat.label_x, "frame");
    sprintf(xpm_mat.label_y, "Fragment");
    sprintf(xpm_mat.legend, "%s", fragment_set->setname);

    snew(xpm_mat.matrix, fastaSequenceSet->nSequences);
    for(i = 0; i < fastaSequenceSet->nSequences; ++i)
        snew(xpm_mat.matrix[i], sequenceLength);

    c.c2 = 0;
    for(i = 0; i < fastaSequenceSet->nSequences; ++i)
        for(j = 0; j < sequenceLength; ++j) {
            c.c1 = fastaSequenceSet->sequence[i].res[j];
            xpm_mat.matrix[i][j] = max(0, searchcmap(xpm_mat.nmap, xpm_mat.map, c));
        }

    write_xpm_m(xpmout,xpm_mat);

    free_xpm_matrix(&xpm_mat, fastaSequenceSet->nSequences, fragment_set->nFragment);
}

/*____________________________________________________________________________*/
/** calculate profile */
void calculate_profile(const char *outFilename, int sequenceLength, SeqSet *fastaSequenceSet, FragmentSet *fragment_set, gmx_bool xpmoutput) {

    int i,j = 0; /* index */
    char *iString; /* string placeholders */
    Set iCodeSet; /* alphabet code set for position i */
    FILE *outprofile = 0; /* file handle for entropy data */

    /* initialize codeSets */
    iCodeSet.nElements = 0;
    iCodeSet.element = safe_malloc(sizeof(Element));

    initialize_code_from_string(fragment_set->codeOrder, &iCodeSet);

    if (my_rank == 0) {
        outprofile = safe_open(outFilename, "w");
        fprintf(outprofile, "         ");
        for(j = 0; j < fragment_set->nFragment; ++j) {
            fprintf(outprofile, " %8c", fragment_set->codeOrder[j]);
        }
        fprintf(outprofile, "\n");
    }

    for (i = 0; i < sequenceLength; ++ i) {
        if (my_rank == 0)
            fprintf(outprofile, " %8d", (i + 1));
        iString = get_string_from_column(fastaSequenceSet, i);
        record_probabilities_from_string(iString, &iCodeSet);

        if (my_rank == 0) {
            for(j = 0; j < iCodeSet.nElements; ++j) {
                fprintf(outprofile, " %8.3f", iCodeSet.element[j].prob);
            }
            fprintf(outprofile, "\n");
        }

        free(iString);
        reset_probabilities(&iCodeSet);

    }

    free(iCodeSet.element);
    fclose(outprofile);

}

/*____________________________________________________________________________*/
/** calculate position entropy */
void entropy_per_position(const char *outFilename, int sequenceLength, SeqSet *fastaSequenceSet, const output_env_t oenv) {

    int i = 0; /* index */
    char *iString; /* string placeholders */
    Set iCodeSet; /* alphabet code set for position i */
    float H = 0.0; /* Shannon's entropy */
    FILE *outEntropy = 0; /* file handle for entropy data */

    if (my_rank == 0)
    outEntropy = xvgropen(outFilename, "Entropy", "Fragment", "H / bit", oenv);
    for (i = 0; i < sequenceLength; ++ i) {
        iString = get_string_from_column(fastaSequenceSet, i);
        extract_code_from_string(iString, &iCodeSet);
        H = Shannon(&iCodeSet);
        if (my_rank == 0)
            fprintf(outEntropy, "%5d %8.3f\n", (i + 1), H);
        free(iCodeSet.element);
        free(iString);
    }
    if (my_rank == 0)
        fclose(outEntropy);
}

/*____________________________________________________________________________*/
/** build data */
void build_data(Data *data, double *xvgvector, char *column, int nSequences, Set *codeSet, FragmentSet *fragment_set, float *maxValue, float *minValue) {

    int i,j,s;
    int allocated = 1;
    int added = 0;
    int (*fcmp)() = &compare_elements;

    codeSet->element = safe_malloc(allocated * sizeof(Element));
    codeSet->nElements = 0;

    *maxValue = FLT_MIN;
    *minValue = FLT_MAX;

    for(i = 0; i < nSequences; ++i) {
        data->instance[i].value = xvgvector[i];
        if (xvgvector[i] > *maxValue)
            *maxValue = xvgvector[i];
        if (xvgvector[i] < *minValue)
            *minValue = xvgvector[i];

        for(s = 0; s < fragment_set->nFragment; ++s)
            if (fragment_set->codeOrder[s] == column[i])
                data->instance[i].code = s;
        data->instance[i].altCode = 0;

        for(j = 0; j < codeSet->nElements; ++j ) {
            if (data->instance[i].code == codeSet->element[j].code) {
                codeSet->element[j].prob += (1.0 / nSequences);
                added = 1;
            }
        }

        if (added == 0) {
            codeSet->element[codeSet->nElements].code = data->instance[i].code;
            codeSet->element[codeSet->nElements].prob = (1.0 / nSequences);
            codeSet->nElements ++;

            if (codeSet->nElements == allocated) {
                allocated += 1;
                codeSet->element = safe_realloc(codeSet->element, allocated * sizeof(Element));
            }
        }

        added = 0;
    }

    MergeSort(codeSet->element, codeSet->nElements, sizeof(Element), fcmp);

    free(column);
}

/*____________________________________________________________________________*/
/** main */
int main(int argc,char *argv[])
{
    /*________________________________________________________________________*/
    /* static variables */
    const char *desc[] = {
        "[TT]g_sa_analyze[tt] performs advanced analyses on the alignment of structural",
        "strings ([TT]-sa[tt]) produced by [TT]g_sa_encode[tt] (Pandini et al., FASEB J. [BB]26[bb], 868",
        "(2012)).[PAR]",

        "Option [TT]-MImatrix[tt] calculates the Mutual Information (MI) matrix between pairs",
        "of columns (positions) in the alignment. The MI values are written to a ASCII",
        "[TT].out[tt] file ([TT]-MImat[tt]). The normalized MI matrix ([TT]-nMImat[tt]), the joint entropy",
        "matrix ([TT]-jHmat[tt]) and the expected finite size error matrix ([TT]-eeMImat[tt]) are also",
        "calculated.[PAR]",

        "The calculation of the statistically significance for the MI values can be requested",
        "with [TT]-nSample[tt] > 0. A random background distribution of [TT]nSample[tt]",
        "samples is generated by shuffling the letters in each column of the",
        "alignment. MI values averaged over the randomized samples ([TT]-meanMImat[tt]), the",
        "associated standard deviations ([TT]-stdMImat[tt]), Z-scores ([TT]-ZMImat[tt]) and [TT]p-values[tt]",
        "([TT]-pvalueMImat[tt]) are evaluated and written to ASCII [TT].out[tt] files.[PAR]",

        "Option [TT]-trmat[tt] writes the transition probability matrix between pairs of SA",
        "letters to an ASCII [TT].out[tt] file ([TT]-trans[tt]).[PAR]",

        "A functional analysis can be requested by providing an [TT].xvg[tt] file with a",
        "time-dependent function-related property ([TT]-fvalue[tt]). The MI values between",
        "selected columns (positions) of the alignment and a discretized version of",
        "the functional property are calculated. The positions in the alignment can be",
        "selected providing a [TT].ndx[tt] file ([TT]-n[tt]). The column containing the functional",
        "property in the [TT].xvg[tt] file can be specified with the [TT]-xvgcol[tt] flag. A summary",
        "of MI statistics is written for each selected position to a separate [TT].log[tt]",
        "file ([TT]-MIlog[tt]). For visualization purposes, the discretized form of",
        "the functional property and a numeric version of the structural string are",
        "written for each selected position to a separate [TT].xvg[tt] file ([TT]-MIxvg[tt]) and to a",
        "separate [TT].out[tt] file ([TT]-MIout[tt]).[PAR]",

        "The same basic statistical analyses performed by [TT]g_sa_encode[tt] can also be",
        "requested:[PAR]",

        "Option [TT]-entropy[tt] writes the Shannon entropy of each position to an [TT].xvg[tt] file",
        "([TT]-Hlf[tt] for local fit or [TT]-Hgf[tt] for global fit).[PAR]",

        "Option [TT]-profile[tt] writes the relative frequency of each letter at each position",
        "to an ASCII [TT].dat[tt] file ([TT]-prolf[tt] for local fit or [TT]-progf[tt] for global fit).[PAR]",

        "Option [TT]-xpm[tt] writes the alignment in [TT].xpm[tt] format: this can be converted to",
        "postscript with [TT]xpm2ps[tt] for an immediate visualization of the time evolution",
        "of the encoding for each position. The name of the [TT].xpm[tt] file can be specified",
        "with the [TT]-xpm[tt] flag."
    };

    static char* userCodeString =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    static gmx_bool entropy = FALSE;
    static gmx_bool profile= FALSE;
    static gmx_bool transmat = FALSE;
    static gmx_bool MImatrix = FALSE;
    static gmx_bool xpmoutput= FALSE;
    static gmx_bool verbose = FALSE;
    static gmx_bool bIndex = FALSE;
    int alphabet=0;
    static char *alphabetname[] = {NULL, "M32K25", NULL};
    enum {a_null, a_M32K25, a_nr};
    Str *fragment_str; /** structure array for alphabet fragments */
    FragmentSet *fragment_set; /** data structure for fragment set */
    FILE *fastaFile = 0; /* fasta local encoding  filehandle */
    FILE *inputFile = 0; /* input local encoding  filehandle */
    FILE *outMatFile = 0; /* ouput trans mat filehandle */
    FILE *xpmout = 0; /* xpm filehandle */
    const char *valueFileName = 0; /* value filename */
    SeqSet fastaSequenceSet; /* fasta sequence set */
    SeqSet inputSequenceSet; /* input sequence set */
    Set iCodeSet; /* alphabet code set for position i */
    ProbMatrix transition_matrix; /* transition matrix */
    int sequenceLength = 0; /* length of encoded sequence */
    char *iString; /* string placeholders */
    SA_colour *SA_colour_set; /* placeholder for colour set  */
    int nColourSets; /* n colour set available */
    double **xvgdata; /* double matrix for xvg values */
    real *outxvgdata[4]; /* double matrix for xvg values output */
    int nx, ny; /* number of x and y */
    int xvgcolumn = 1; /* xvg column */
    float binwidth = 0.01; /* bin width */
    float *breaks = 0; /* break array */
    int *breaks_idx = 0; /* break index array */
    Collection bins; /* bin data structure */
    int nBins; /* number of bins */
    gmx_bool notAlloc = FALSE; /* flag for allocation */
    Data data; /* values */
    Set codeSet; /* code set */
    Set altCodeSet; /* alt code set */
    float maxValue, minValue; /* max and min values */
    char basenameMI[64]; /* MI basename */
    char MIxvgFilename[64]; /* MI xvg filename */
    char MIoutFilename[64]; /* MI out filename */
    char MIlogFilename[64]; /* MI log filename */
    const char *MIxvglegend[] = {"value", "SA", "partition"}; /* MI xvg legend */
    ProbMatrix probMat; /* probability matrix */
    float mutualInformation; /* MI */
    float jointEntropy; /* joint entropy */
    FILE *MIlog; /* MI log filehandle */
    FILE *MIout; /* MI out filehandle */
    int nSurrogates = 0; /* n random samples for MI Matrix significance */
    Set jCodeSet; /* alphabet code set for position j */
    float **MIMat = 0; /* Mutual Information matrix */
    float **JentropyMat = 0; /* Joint entropy matrix */
    float **eeMIMat = 0; /* estimated Mutual Information error matrix */
    float eeMI; /* estimated Mutual Information error */
    float meanMI = 0.0; /* mean random Mutual Information */
    float stdMI = 0.0; /* std random Mutual Information */
    float pValueMI = 0.0; /*  p-value of Mutual Information */
    float **meanMIMat = 0; /* mean random Mutual Information matrix */
    float **stdMIMat = 0; /* std random Mutual Information matrix */
    float **pValueMIMat = 0; /* p-value Mutual Information matrix */
    time_t now; /* actual time */
    gsl_rng *rndGenerator = 0; /* random number generator */
    const gsl_rng_type *rndGeneratorType; /* random number generator type */
    FILE *MIFile = 0; /* MI matrix file handle */
    FILE *eeMIFile = 0; /* error matrix file handle */
    FILE *jHFile = 0; /* joint entropy matrix file handle */
    FILE *NormMIFile = 0; /* normalised MI matrix file handle */
    FILE *ZMIFile = 0; /* Z-score MI matrix file handle */
    FILE *MeanMIFile = 0; /* mean MI matrix file handle */
    FILE *StdMIFile = 0; /* std deviation MI matrix file handle */
    FILE *pValueMIFile = 0; /* p-value MI matrix file handle */
    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* initialize MPI routines */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    /*________________________________________________________________________*/
    /* Extra arguments */

    t_pargs pa[] = {
        {   "-alphabet", FALSE, etENUM, {&alphabetname},
            "Structural Alphabet"
        },
        {   "-entropy", FALSE, etBOOL, {&entropy},
            "Calculate entropy per position"
        },
        {   "-profile", FALSE, etBOOL, {&profile},
            "Calculate frequency profile"
        },
        {   "-trmat", FALSE, etBOOL, {&transmat},
            "Calculate transition matrix"
        },
        {   "-MImatrix", FALSE, etBOOL, {&MImatrix},
            "Calculate Mutual Information matrix"
        },
        {   "-nSample", FALSE, etINT, {&nSurrogates},
            "Number of random samples for significance estimation"
        },
        {   "-xpm", FALSE, etBOOL, {&xpmoutput},
            "Output xpm matrix"
        },
        {   "-xvgcol", FALSE, etINT, {&xvgcolumn},
            "xvg column"
        },
        {   "-bw", FALSE, etREAL, {&binwidth},
            "bin width for partition"
        },
        {   "-verbose", FALSE, etBOOL, {&verbose},
            "Verbose output"
        }
    };

    /*________________________________________________________________________*/
    /* gromacs variables */

    t_topology top;
    t_atoms    *atoms=NULL;
    int        ePBC;
    char       title[STRLEN];
    t_trxframe fr;
    rvec       *xtop;
    matrix     box;
    int        status;
    int        flags = TRX_READ_X;
    int        i,j;
    char      *tgrpname = 0;
    char      *vgrpname = 0;
    int       isize = 0;
    atom_id   *index = 0;
    output_env_t oenv;

    /*________________________________________________________________________*/
    /* gromacs file variables */

    t_filenm fnm[] = {
        { efOUT, "-sa", "lf_str.out", ffREAD},     /* sa input */
        { efNDX, NULL, "pos_index",  ffOPTRD},     /* index file */
        { efXVG, "-fvalue", "funval", ffOPTRD},    /* xvg input values */
        { efXVG, "-H", "lf_entropy", ffOPTWR},     /* entropy */
        { efDAT, "-pro", "lf_prof", ffOPTWR},      /* profile */
        { efOUT, "-trans", "transmat", ffOPTWR},   /* trans mat output */
        { efXPM, "-xpmout", "lf", ffOPTWR},        /* xpm matrix */
        { efXVG, "-MIxvg", "lf_MI", ffOPTWR},      /* MI xvg output*/
        { efLOG, "-MIlog", "lf_MI", ffOPTWR},      /* MI log output*/
        { efOUT, "-MIout", "lf_MI", ffOPTWR},      /* MI out output*/
        { efOUT, "-MImat", "lf_MImat", ffOPTWR},       /* MI matrix output*/
        { efOUT, "-eeMImat", "lf_eeMImat", ffOPTWR},       /* error matrix output*/
        { efOUT, "-jHmat", "lf_jHmat", ffOPTWR},       /* joint entropy output*/
        { efOUT, "-nMImat", "lf_nMImat", ffOPTWR},       /* normalised MI output*/
        { efOUT, "-ZMImat", "lf_ZMImat", ffOPTWR},       /* Z-score MI matrix file handle */
        { efOUT, "-meanMImat", "lf_meanMImat", ffOPTWR},       /* mean MI matrix file handle */
        { efOUT, "-stdMImat", "lf_stdMImat", ffOPTWR},       /* std deviation MI matrix file handle */
        { efOUT, "-pvalueMImat", "lf_pvalueMImat", ffOPTWR}       /* p-value MI matrix file handle */
    };

#define NFILE asize(fnm)

    /*________________________________________________________________________*/
    /* Print (C) */
    CopyRight(stderr,argv[0]);

    /*________________________________________________________________________*/
    /* Parse args */
    parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm,
                      asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);

    /*________________________________________________________________________*/
    /* check alphabet */
    alphabet=1;
    while (alphabet < a_nr &&
            strcasecmp(alphabetname[0], alphabetname[alphabet])!= 0)
        alphabet++;
    if (alphabet == a_nr)
        gmx_fatal(FARGS,"Invalid alphabet");

    if (my_rank == 0)
        fprintf(stderr,"Using %s alphabet for encoding\n", alphabetname[0]);

    /*________________________________________________________________________*/
    /** Load the fragment set data into fragment_set and fragment_str */
    fragment_set = safe_malloc(sizeof(FragmentSet));
    fragment_str = load_fragment_data(0, userCodeString, alphabetname[0],
                                      fragment_set, 0);

    /*________________________________________________________________________*/
    /** open input file */
    inputFile = opt2FILE("-sa", asize(fnm), fnm, "r");
    /** read input file */
    read_inputFileEnsemble(inputFile, &inputSequenceSet);
    /** close input file */
    fclose(inputFile);

    /*________________________________________________________________________*/
    /* initialize code set */
    iCodeSet.nElements = 0;
    iCodeSet.element = safe_malloc(sizeof(Element));

    initialize_code_from_string(fragment_set->codeOrder, &iCodeSet);

    /* sequence length */
    sequenceLength = inputSequenceSet.sequence[0].length;

    /*________________________________________________________________________*/
    /* calculate per position entropy */
	if (my_rank == 0) fprintf(stdout, "Entropy per position\n");
    if (entropy)
        entropy_per_position(opt2fn("-H", asize(fnm), fnm), sequenceLength, &inputSequenceSet, oenv);

    /*________________________________________________________________________*/
    /* calculate profile */
	if (my_rank == 0) fprintf(stdout, "Profile\n");
    if (profile)
        calculate_profile(opt2fn("-pro", asize(fnm), fnm), sequenceLength, &inputSequenceSet, fragment_set, xpmoutput);

    /*________________________________________________________________________*/
    /* calculate transition matrix */
	if (my_rank == 0) fprintf(stdout, "Transition matrix\n");
    if (transmat) {

        /*________________________________________________________________________*/
        /* read index file if provided */
        bIndex = ftp2bSet(efNDX,NFILE,fnm);
        if (bIndex) {
            if (my_rank == 0)
                fprintf(stderr, "\nSelect group for transition matrix calculation:\n");
            rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&tgrpname);
        }

        /*________________________________________________________________________*/
        /** calculate transition matrix */
        /* initialize the matrix */
        initialize_probability_matrix(&transition_matrix, &iCodeSet, &iCodeSet);

        /* sequence length */
        sequenceLength = inputSequenceSet.sequence[0].length;

        calculate_transition_frequencies(&inputSequenceSet, &transition_matrix, isize, index);

        if (my_rank == 0) {
            outMatFile = opt2FILE("-trans", asize(fnm), fnm, "w");
            fprintf(outMatFile, "          ");
            for(j = 0; j < fragment_set->nFragment; ++j)
                fprintf(outMatFile, "%10c", transition_matrix.altCodeSet->element[j].code);
            fprintf(outMatFile, "\n");
            for(i = 0; i < fragment_set->nFragment; ++i) {
                fprintf(outMatFile, "%10c", transition_matrix.codeSet->element[i].code);
                for(j = 0; j < fragment_set->nFragment; ++j)
                    fprintf(outMatFile, "%10.5f", transition_matrix.prob[i][j]);
                fprintf(outMatFile, "\n");
            }
            fclose(outMatFile);
        }

        if (bIndex) {
            sfree(index);
            sfree(tgrpname);
        }
    }

    /*________________________________________________________________________*/
    /* print xpm map of encoding */
    SA_colour_set = load_colour_set(&nColourSets);

    if (xpmoutput) {
        xpmout = opt2FILE("-xpmout", asize(fnm), fnm, "w");
        save_SA_xpm_matrix(xpmout, &inputSequenceSet, sequenceLength, fragment_set, SA_colour_set, nColourSets);
        fclose(xpmout);
    }

    /*________________________________________________________________________*/
    /* read value file if provided */
    valueFileName = opt2fn_null("-fvalue", NFILE, fnm);

    if (valueFileName != 0) {

        if (binwidth <= 0) {
            if (my_rank == 0)
                fprintf(stderr, "Warning: bin width should be positive.\n");
            notAlloc = TRUE;
        } else {
            nx = read_xvg(valueFileName, &xvgdata, &ny);

            if (nx != inputSequenceSet.nSequences) {
                if (my_rank == 0)
                    fprintf(stderr, "Warning: number of value differs from number of frames.\n");
                notAlloc = TRUE;
            } else {
                if (xvgcolumn > (ny - 1)) {
                    if (my_rank == 0)
                        fprintf(stderr, "Warning: xvgcol number larger than number of columns.\n");
                    notAlloc = TRUE;
                } else {

                    /* allocate memory for output vectors */
                    outxvgdata[0] = safe_malloc(sizeof(real) * nx);
                    outxvgdata[1] = safe_malloc(sizeof(real) * nx);
                    outxvgdata[2] = safe_malloc(sizeof(real) * nx);
                    outxvgdata[3] = safe_malloc(sizeof(real) * nx);

                    /* read index file if provided */
                    bIndex = ftp2bSet(efNDX,NFILE,fnm);
                    if (bIndex) {
                        if (my_rank == 0)
                            fprintf(stderr, "\nSelect group for value partition:\n");
                        rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&vgrpname);
                    }

                    /** allocate memory for bins */
                    for(i = 0; i * binwidth < 1.00; i++)
                        continue;
                    nBins = i + 1;
                    bins.object = safe_malloc(sizeof(Object) * (nBins + 1));
                    bins.nObjects = nBins;

                    for(i = 0; i < nBins; i++)
                        bins.object[i + 1].value = i * binwidth;

                    /** link data to bins */
                    bins.data = 0;

                    /** allocate memory for breaks data */
                    breaks = safe_malloc(fragment_set->nFragment * sizeof(float));
                    breaks_idx = safe_malloc(fragment_set->nFragment * sizeof(int));

                    if (isize == 0) {
                        if (my_rank == 0) {
                            fprintf(stderr, "\nWarning: a value file has been provided but it will not be used\n");
                            fprintf(stderr, "         because no index file has been provided to select the\n");
                            fprintf(stderr, "         group of fragments for the analysis.\n");
                        }
                    }
                    for(i = 0; i < isize; ++ i) {

                        if (my_rank == 0)
                            fprintf(stderr, "Now processing... position %d\n", index[i] + 1);

                        strncpy(basenameMI, opt2fn("-MIxvg", asize(fnm), fnm), strlen(opt2fn("-MIxvg", asize(fnm), fnm)) - 4);
                        basenameMI[strlen(opt2fn("-MIxvg", asize(fnm), fnm)) - 4] = '\0';

                        /** allocate memory for data */
                        data.instance = malloc(sizeof(Instance) * nx);
                        data.nInstances = nx;

                        /** build data */
                        build_data(&data, xvgdata[xvgcolumn],
                                   get_string_from_column(&inputSequenceSet, index[i]),
                                   inputSequenceSet.nSequences, &codeSet,
                                   fragment_set, &maxValue, &minValue);
                        for(j = 0; j < data.nInstances; ++j) {
                            outxvgdata[0][j] = xvgdata[0][j];
                            outxvgdata[1][j] = data.instance[j].value;
                            outxvgdata[2][j] = (real) data.instance[j].code;
                            data.instance[j].value = (data.instance[j].value - minValue) / (maxValue - minValue);
                        }
                        data.codeSet = &codeSet;
                        bins.data = &data;

                        /** partition data */
                        partition(&bins, nBins, fragment_set->nFragment, breaks, breaks_idx, 1, 2);

                        /** update alternate code */
                        getAltCode(breaks, fragment_set->nFragment, &data, &altCodeSet);

                        for(j = 0; j < data.nInstances; ++j)
                            outxvgdata[3][j] = (real) bins.data->instance[j].altCode;

                        /** write output xvg table */
                        if (my_rank == 0) {
                            if (isize > 1) {
                                sprintf(MIxvgFilename, "%s.%04d.xvg", basenameMI, index[i] + 1);
                                write_xvg(MIxvgFilename, "MI Partition", nx, 4, outxvgdata, MIxvglegend, oenv);
                            } else {
                                write_xvg(opt2fn("-MIxvg", asize(fnm), fnm), "MI Partition", nx, 4, outxvgdata, MIxvglegend, oenv);
                            }
                        }

                        /** write output log table */
                        if (my_rank == 0) {
                            if (isize > 1) {
                                sprintf(MIoutFilename, "%s.%04d.out", basenameMI, index[i] + 1);
                                MIout = safe_open(MIoutFilename, "w");
                            } else {
                                MIout = safe_open(opt2fn("-MIout", asize(fnm), fnm), "w");
                            }
                            for(j = 0; j < data.nInstances; ++j)
                                fprintf(MIout, "%10.3f %8.3f %8.3f %c %2d\n",
                                        xvgdata[0][j], xvgdata[xvgcolumn][j], data.instance[j].value,
                                        fragment_set->codeOrder[data.instance[j].code], data.instance[j].altCode);
                            fclose(MIout);
                        }

                        /** initialize and populate probability matrix */
                        initialize_probability_matrix(&probMat, &codeSet, &altCodeSet);
                        populate_probability_matrix(&probMat, &data);
                        mutualInformation = mutual_information(&probMat, &codeSet, &altCodeSet);
                        jointEntropy = joint_entropy(&probMat, &codeSet, &altCodeSet);

                        /** print entropy and mutual information */
                        if (my_rank == 0) {
                            if (isize > 1) {
                                sprintf(MIlogFilename, "%s.%04d.log", basenameMI, index[i] + 1);
                                MIlog = safe_open(MIlogFilename, "w");
                            } else {
                                MIlog = safe_open(opt2fn("-MIlog", asize(fnm), fnm), "w");
                            }

                            fprintf(MIlog, "      Code Entropy: %10.5f\n", Shannon(&codeSet));
                            fprintf(MIlog, "   altCode Entropy: %10.5f\n", Shannon(&altCodeSet));
                            fprintf(MIlog, "Mutual Information: %10.5f\n", mutualInformation);
                            fprintf(MIlog, "     joint Entropy: %10.5f\n", jointEntropy);
                            /*
                            fprintf(MIlog, "        quantity D: %10.5f\n", 1 - (mutualInformation / jointEntropy));
                            if (Shannon(&codeSet) > Shannon(&altCodeSet))
                            	fprintf(MIlog, "       quantity D': %10.5f\n", 1 - (mutualInformation / Shannon(&codeSet)));
                            else
                            	fprintf(MIlog, "       quantity D': %10.5f\n", 1 - (mutualInformation / Shannon(&altCodeSet)));
                            */

                            fclose(MIlog);

                            /** free memory from data */
                            free(data.instance);
                            free_probability_matrix(&probMat);
                        }
                    }

                    /** free index an grpname memory */
                    if (bIndex) {
                        sfree(index);
                        sfree(vgrpname);
                    }
                }
            }
        }
    } else {
        bIndex = ftp2bSet(efNDX,NFILE,fnm);
        if ((!transmat) && bIndex) {
            if (my_rank == 0) {
                fprintf(stderr, "\nWarning: an index file has been provided but it will not be used.\n");
                fprintf(stderr, "         This type of file is only used with the option -trmat or\n");
                fprintf(stderr, "         the option -fvalue.\n");
            }
        }
    }

    if (MImatrix) {
        /*________________________________________________________________________*/
        /* initialize code set */
        jCodeSet.nElements = 0;
        jCodeSet.element = safe_malloc(sizeof(Element));

        initialize_code_from_string(fragment_set->codeOrder, &jCodeSet);

        /* reset probabilities */
        for(i = 0; i < iCodeSet.nElements; ++i) {
            iCodeSet.element[i].prob = 0.0;
            jCodeSet.element[i].prob = 0.0;
        }

        /*________________________________________________________________________*/
        /** allocate and initialize Mutual Information matrix */
        MIMat = alloc_float_matrix(MIMat, sequenceLength, sequenceLength);
        initialise_float_matrix(MIMat, sequenceLength, sequenceLength, 0.0);
        eeMIMat = alloc_float_matrix(eeMIMat, sequenceLength, sequenceLength);
        initialise_float_matrix(eeMIMat, sequenceLength, sequenceLength, 0.0);
        JentropyMat = alloc_float_matrix(JentropyMat, sequenceLength, sequenceLength);
        initialise_float_matrix(JentropyMat, sequenceLength, sequenceLength, 0.0);

        /*________________________________________________________________________*/
        /** allocate and initialize probability matrix */
        initialize_probability_matrix(&probMat, &iCodeSet, &jCodeSet);

        /*________________________________________________________________________*/
        /** calculate Mutual Information and Joint entropy */
		int zz;
		int completion_i = 0;
		double completion;

		MPI_Barrier(MPI_COMM_WORLD);

		if (my_rank == 0) fprintf(stdout, "\nMutual Information and Joint Entropy\n");
		fflush(stdout);
        for(i = 0, zz = 0, completion_i = 0; i < sequenceLength; ++i) {
            for(j = i; j < sequenceLength; ++j) {
			   /* print progress */
				++ zz; 
				completion = (long double)zz / ((long double)(sequenceLength*(sequenceLength-1)) / 2) * 100;
				if ((int)completion > completion_i) {
					completion_i = (int)completion;
					if (my_rank == 0) {
						fprintf(stdout, "\t%3d%%\r", completion_i);
						fflush(stdout);
					}
				}

                MIMat[i][j] = column_mutual_information(&inputSequenceSet, &probMat, i, j, &eeMI);
                eeMIMat[i][j] = eeMI;
                if (i != j) {
                    MIMat[j][i] = MIMat[i][j];
                    eeMIMat[j][i] = eeMIMat[i][j];
                }
            }
        }
        for(i = 0; i < sequenceLength; ++i) {
            for(j = i; j < sequenceLength; ++j) {
                JentropyMat[i][j] = column_joint_entropy(&inputSequenceSet, &probMat, i, j);
                if (i != j)
                    JentropyMat[j][i] = JentropyMat[i][j];
            }
        }

        /*________________________________________________________________________*/
        /** output Mutual Information matrices */
        /** output MI */
        if (my_rank == 0) {
            MIFile = safe_open(opt2fn("-MImat", asize(fnm), fnm), "w");
            for(i = 0; i < sequenceLength; ++i) {
                for(j = 0; j < sequenceLength; ++j) {
                    fprintf(MIFile, "%8.3f", MIMat[i][j]);
                }
                fprintf(MIFile, "\n");
            }
            fclose(MIFile);
        }

        /** output error MI */
        if (my_rank == 0) {
            eeMIFile = safe_open(opt2fn("-eeMImat", asize(fnm), fnm), "w");
            for(i = 0; i < sequenceLength; ++i) {
                for(j = 0; j < sequenceLength; ++j) {
                    fprintf(eeMIFile, "%8.3f", eeMIMat[i][j]);
                }
                fprintf(eeMIFile, "\n");
            }
            fclose(eeMIFile);
        }

        /** output joint H */
        if (my_rank == 0) {
            jHFile = safe_open(opt2fn("-jHmat", asize(fnm), fnm), "w");
            /** write output joint H file */
            for(i = 0; i < sequenceLength; ++i) {
                for(j = 0; j < sequenceLength; ++j) {
                    fprintf(jHFile, "%8.3f", JentropyMat[i][j]);
                }
                fprintf(jHFile, "\n");
            }
            /** close output joint H file */
            fclose(jHFile);
        }

        /** output normalized MI */
        if (my_rank == 0) {
            NormMIFile = safe_open(opt2fn("-nMImat", asize(fnm), fnm), "w");
            for(i = 0; i < sequenceLength; ++i) {
                for(j = 0; j < sequenceLength; ++j) {
                    fprintf(NormMIFile, "%8.3f", (MIMat[i][j] / JentropyMat[i][j]));
                }
                fprintf(NormMIFile, "\n");
            }
            fclose(NormMIFile);
        }

        /*________________________________________________________________________*/
        /** perform significance analysis */
		if (my_rank == 0) fprintf(stdout, "\nSignificance analysis\n");
		fflush(stdout);
        if (nSurrogates > 0) {
            /*________________________________________________________________________*/
            /** initialize random number generator */
            /** get time */
            now = time(0);
            /** initialize random number generator environment */
            gsl_rng_env_setup();
            rndGeneratorType = gsl_rng_default;
            rndGenerator = gsl_rng_alloc(rndGeneratorType);
            /** seed random number generator */
            gsl_rng_set(rndGenerator, (now + my_rank));

            /*________________________________________________________________________*/
            /** allocate and initialize matrices for significance analysis */
            meanMIMat = alloc_float_matrix(meanMIMat, sequenceLength, sequenceLength);
            initialise_float_matrix(meanMIMat, sequenceLength, sequenceLength, 0.0);
            stdMIMat = alloc_float_matrix(stdMIMat, sequenceLength, sequenceLength);
            initialise_float_matrix(stdMIMat, sequenceLength, sequenceLength, 0.0);
            pValueMIMat = alloc_float_matrix(pValueMIMat, sequenceLength, sequenceLength);
            initialise_float_matrix(pValueMIMat, sequenceLength, sequenceLength, 0.0);

            /*____________________________________________________________________*/
            /** calculate Mutual Information surrogate matrices */
			MPI_Barrier(MPI_COMM_WORLD);

            for(i = 0, zz = 0, completion_i = 0; i < sequenceLength; ++i) {
                for(j = i; j < sequenceLength; ++j) {
				   /* print progress */
					++ zz; 
					completion = (long double)zz * nodes / ((long double)(sequenceLength*(sequenceLength-1)) / 2) * 100;
					if ((int)completion > completion_i) {
						completion_i = (int)completion;
						fprintf(stdout, "\t%3d%%\r", completion_i);
						fflush(stdout);
					}

                    random_column_mutual_information(&meanMI, &stdMI, &pValueMI, MIMat[i][j], &inputSequenceSet, &probMat, i, j, nSurrogates, rndGenerator);
                    meanMIMat[i][j] = meanMI;
                    stdMIMat[i][j] = stdMI;
                    pValueMIMat[i][j] = pValueMI;

                    if (i != j) {
                        meanMIMat[j][i] = meanMIMat[i][j];
                        stdMIMat[j][i] = stdMIMat[i][j];
                        pValueMIMat[j][i] = pValueMIMat[i][j];
                    }
                }
            }

            /*____________________________________________________________________*/
            /** output significance data */
            /** open output MI files */
            if (my_rank == 0) {
                ZMIFile = safe_open(opt2fn("-ZMImat", asize(fnm), fnm), "w");
                MeanMIFile = safe_open(opt2fn("-meanMImat", asize(fnm), fnm), "w");
                StdMIFile = safe_open(opt2fn("-stdMImat", asize(fnm), fnm), "w");
                pValueMIFile = safe_open(opt2fn("-pvalueMImat", asize(fnm), fnm), "w");

                /** write output MI files */
                for(i = 0; i < sequenceLength; ++i) {
                    for(j = 0; j < sequenceLength; ++j) {
                        /* set significance of I(X,X) to 0.0 and p-value 0.50 */
                        if (i == j) {
                            fprintf(ZMIFile, "%10.3f", 0.0);
                            fprintf(MeanMIFile, "      NA");
                            fprintf(StdMIFile, "      NA");
                            fprintf(pValueMIFile, "%16.5f", 0.50);
                        } else {
                            fprintf(ZMIFile, "%10.3f", ((MIMat[i][j] - meanMIMat[i][j]) / stdMIMat[i][j]) );
                            fprintf(MeanMIFile, "%8.3f", meanMIMat[i][j]);
                            fprintf(StdMIFile, "%8.3f", stdMIMat[i][j]);
                            fprintf(pValueMIFile, "%16.5f", pValueMIMat[i][j]);
                        }
                    }
                    fprintf(ZMIFile, "\n");
                    fprintf(MeanMIFile, "\n");
                    fprintf(StdMIFile, "\n");
                    fprintf(pValueMIFile, "\n");
                }

                /** close output MI files */
                fclose(ZMIFile);
                fclose(MeanMIFile);
                fclose(StdMIFile);
                fclose(pValueMIFile);
            }

        }
    }

    /*________________________________________________________________________*/
    /** close file handles */

    /*________________________________________________________________________*/
    /* free memory */

    /* structure */
    for (i = 0; i < fragment_set->nFragment; ++ i)
        free(fragment_str[i].atom);
    free(fragment_str);

    /* fragment set */
    free_float_matrix3D(fragment_set->coord_values, fragment_set->nFragment, fragment_set->lFragment);
    free(fragment_set);

    /* colours */
    free_colour_set(SA_colour_set);

    /* code set */
    free(iCodeSet.element);
    if (MImatrix)
        free(jCodeSet.element);

    /* code set */
    if (transmat)
        free_probability_matrix(&transition_matrix);

    /* input sequence set */
    for(i = 0; i < inputSequenceSet.nSequences; ++i) {
        free(inputSequenceSet.sequence[i].res);
        free(inputSequenceSet.sequence[i].name);
    }
    free(inputSequenceSet.sequence);

    /** partition data structures */
    if ((valueFileName != 0) & (!notAlloc)) {
        free(codeSet.element);
        free(altCodeSet.element);
        /** collection */
        free(bins.object);
        /** breaks */
        free(breaks);
        free(breaks_idx);

        /* xvg data */
        for(i = 0; i < ny; ++i)
            sfree(xvgdata[i]);
        sfree(xvgdata);

        /* xvg output data */
        free(outxvgdata[0]);
        free(outxvgdata[1]);
        free(outxvgdata[2]);
        free(outxvgdata[3]);
    }

    if (MImatrix) {
        /** free Mutual Information matrix */
        free_float_matrix(MIMat, sequenceLength);
        free_float_matrix(eeMIMat, sequenceLength);
        free_float_matrix(JentropyMat, sequenceLength);

        /** free probability matrix */
        free_probability_matrix(&probMat);

    }

    if (MImatrix & (nSurrogates > 0)) {
        /** free  matrices for significance analysis */
        free_float_matrix(meanMIMat, sequenceLength);
        free_float_matrix(stdMIMat, sequenceLength);
        free_float_matrix(pValueMIMat, sequenceLength);

        /** free random number generator */
        gsl_rng_free(rndGenerator);
    }

    /*________________________________________________________________________*/
    /* MPI */
#ifdef MPI
    /* stop MPI processes */
    MPI_Finalize();
#endif

    /*________________________________________________________________________*/
    return 0;
}

