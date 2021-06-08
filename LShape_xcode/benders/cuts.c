/*
 * cuts.c
 *
 *  Created on: Sep 21, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "benders.h"

extern configType config;

#if defined(SAVE_DUALS)
extern dualsType *duals;
#endif

int formOptCut(oneProblem *orig, stocType *stoc, probType *prob, cellType *cell, vector Xvect, BOOL isIncumb) {
	oneCut 	*cut;
	int    	cutIdx, obs,  t, i, idx;
    char *str_obs= malloc(2);
    char *u= malloc(6);

	/* Only a fraction (at least one) of subproblems are solved in any iteration. */
    cut = newCut(prob->num->prevCols, cell->k);

    for (obs=0; obs < 51; obs++){
        sprintf(str_obs, "%d", obs);
        u[0] = '\0';
        strcat(u, "u(");
        strcat(u, str_obs);
        strcat(u, ")");
        for (t=0; t < 7; t++){
            for (i=0; i < 3; i++){
                for (idx = 1; idx <= prob->num->prevCols; idx++)
                    if (!strcmp(orig->cname[idx-1], u))
                        break;
                if (stoc->vals[stoc->groupBeg[obs]+t][i] > Xvect[idx]){
                    cut->alpha += 1.0 / 3 * stoc->vals[stoc->groupBeg[obs]+t][i];
                    cut->beta[idx] += 1.0 / 3;
                    }
                }
            }
        }
    
	cut->beta[0] = 1.0;
	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, prob->num->prevCols);
	else
		cut->alphaIncumb = cut->alpha;
    
#if defined(STOCH_CHECK)
        printf("Objective estimate computed as cut height = %lf\n", cut->alpha - vXv(cut->beta, Xvect, NULL, prob->num->prevCols));
#endif

	/* (c) add cut to the master problem  */
	if ( (cutIdx = addCut2Master(cell, cell->cuts, cut, prob->num->prevCols)) < 0 ) {
		errMsg("algorithm", "formSDCut", "failed to add the new cut to master problem", 0);
		goto TERMINATE;
	}

    mem_free(u);
    mem_free(str_obs);
	return cutIdx;
	TERMINATE: 	return -1;
}//END formCut()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, vector xk, int betaLen) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], xk, betaLen);
		if (Sm < ht)
			Sm = ht;
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X. */
double cutHeight(oneCut *cut, vector xk, int betaLen) {
	double height;

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	return height;
}//END cutHeight()

/* This function will remove the oldest cut whose corresponding dual variable is zero (thus, a cut which was slack in last solution). */
int reduceCuts(oneProblem *master, cutsType *cuts, vector vectX, vector piM, int betaLen, int *iCutIdx) {
	double 	height, minHeight;
	int 	oldestCut, idx;

	oldestCut = cuts->cnt;

	/* identify the oldest loose cut */
	for (idx = 0; idx < cuts->cnt; idx++) {
		if ( idx == (*iCutIdx) || cuts->vals[idx]->rowNum < 0)
			/* avoid dropping incumbent cut and newly added cuts */
			continue;

		if ( cuts->vals[idx]->ck < oldestCut && DBL_ABS(piM[cuts->vals[idx]->rowNum + 1]) <= config.TOLERANCE ) {
			oldestCut = idx;
		}
	}

	/* if the oldest loose cut is the most recently added cut, then the cut with minimium cut height will be dropped */
	if ( oldestCut == cuts->cnt ) {
		minHeight = cutHeight(cuts->vals[0], vectX, betaLen);
		oldestCut = 0;

		for (idx = 1; idx < cuts->cnt; idx++) {
			if (idx == (*iCutIdx))
				continue;

			height = cutHeight(cuts->vals[idx], vectX, betaLen);
			if (height < minHeight) {
				minHeight = height;
				oldestCut = idx;
			}
		}
	}

	/* drop the selected cut and swap the last cut into its place */
	if ( dropCut(master, cuts, oldestCut, iCutIdx) ){
		errMsg("algorithm", "reduceCuts", "failed to drop a cut", 0);
		return -1;
	}

	return oldestCut;
}//END reduceCuts()

/* This function removes a cut from both the cutType structure and the master problem constraint matrix.  In the cuts->vals array, the last
 * cut is swapped into the place of the exiting cut.  In the constraint matrix, the row is deleted, and the row numbers of all constraints
 * below it are decremented. */
int dropCut(oneProblem *master, cutsType *cuts, int cutIdx, int *iCutIdx) {
	int idx, deletedRow;

	deletedRow = cuts->vals[cutIdx]->rowNum;
	/* Get rid of the indexed cut on the solver */
	if (  removeRow(master->lp, deletedRow, deletedRow) ) {
		errMsg("solver", "dropCut", "failed to remove a row from master problem", 0);
		return 1;
	}
	freeOneCut(cuts->vals[cutIdx]);

	/* move the last cut to the deleted cut's position (structure) */
	cuts->vals[cutIdx] = cuts->vals[--cuts->cnt];

	/* if the swapped cut happens to be the incumbent cut, then update its index */
	if ( (*iCutIdx) == cuts->cnt )
		(*iCutIdx) = cutIdx;

	for (idx = 0; idx < cuts->cnt; idx++) {
		if (cuts->vals[idx]->rowNum > deletedRow)
			--cuts->vals[idx]->rowNum;
	}

	/* decrease the number of rows on solver */
	master->mar--;

	return 0;
}//END dropCut()

oneCut *newCut(int numX, int currentIter) {
	oneCut *cut;
    int i;

	cut = (oneCut *) mem_malloc (sizeof(oneCut));
	cut->isIncumb = FALSE; 								/* new cut is by default not an incumbent */
	cut->alphaIncumb = 0.0;
	cut->rowNum = -1;
	cut->ck = currentIter;
	if (!(cut->beta = arr_alloc(numX + 1, double)))
		errMsg("allocation", "new_cut", "beta", 0);
    for (i=0; i < numX + 1; i++){
        cut->beta[i] = 0;
    }
	cut->alpha = 0.0;
	cut->name = (string) arr_alloc(NAMESIZE, char);

	return cut;
}//END newCut

cutsType *newCuts(int maxCuts) {
	cutsType *cuts;

	if (maxCuts == 0)
		return NULL;

	if (!(cuts = (cutsType *) mem_malloc (sizeof(cutsType))))
		errMsg("allocation", "newCuts", "cuts",0);
	if (!(cuts->vals = (oneCut **) arr_alloc (maxCuts, oneCut)))
		errMsg("allocation", "newCuts", "oneCuts",0);
	cuts->cnt = 0;

	return cuts;
}//END newCuts

void freeOneCut(oneCut *cut) {

	if (cut) {
		if (cut->beta) mem_free(cut->beta);
		if (cut->name) mem_free(cut->name);
		mem_free(cut);
	}
}

void freeCutsType(cutsType *cuts, BOOL partial) {
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		freeOneCut(cuts->vals[cnt]);

	if ( partial )
		cuts->cnt = 0;
	else {
		mem_free(cuts->vals);
		mem_free(cuts);
	}
}//END freeCuts
