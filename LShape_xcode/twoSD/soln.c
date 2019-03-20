/*
 * soln.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

/***********************************************************************\
 ** This function determines whether the "stagewise descent property" is
 ** satisified.  If the current approximation of f_k gives a lower difference
 ** between the candidate and incumbent x than the previous approximation
 ** gave, then the incumbent x is updated to the candidate x, and the
 ** reference to the incumbent cut is updated as well.  The function returns
 ** TRUE if the incumbent was updated; FALSE otherwise.
 \***********************************************************************/
int checkImprovement(probType *prob, cellType *cell, int candidCut) {
	double  candidEst;

	/* Calculate height at new candidate x with newest cut included */
	candidEst = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, prob->num->cols, cell->lb);;
	cell->incumbEst = vXvSparse(cell->incumbX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->incumbX, prob->num->cols, cell->lb);

#ifdef ALGO_CHECK
	printf("\nCandidate estimate = %lf, Incumbent estimate = %lf",candidEst, cell->incumbEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ((candidEst - cell->incumbEst) < (config.R1 * cell->gamma)) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(prob, cell, candidEst) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}
		cell->iCutIdx = candidCut;
		cell->incumbChg = FALSE;
		printf("+"); fflush(stdout);
	}
	else {
		/* Update quad_scalar when no incumbent is found. */
		cell->quadScalar = min(config.MAX_QUAD_SCALAR, cell->quadScalar / config.R2);
		cell->normDk_1 = cell->normDk;
	}

	if ( changeQPproximal(cell->master->lp, prob->num->cols, cell->quadScalar) ) {
		errMsg("setup", "checkImprovement", "failed to add the proximal term to QP", 0);
		return 1;
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(probType *prob, cellType *cell, double candidEst) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, 1);
	cell->incumbEst = candidEst;

	/* update the proximal parameter based on estimated improvement */
	if ( cell->normDk > config.TOLERANCE )
		if ( cell->normDk >= config.R3 * cell->normDk_1 ) {
			cell->quadScalar *= config.R2 * config.R3 * cell->normDk_1/ cell->normDk;
			cell->quadScalar  = min(config.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = max(config.MIN_QUAD_SCALAR, cell->quadScalar);
		}

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( constructQP(prob, cell, cell->incumbX, cell->quadScalar) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->iCutUpdt = cell->k;
	cell->incumbChg = TRUE;

	/* keep the two norm of solution*/
	cell->normDk_1 = cell->normDk;
	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = FALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	cell->gamma = 0.0;

	return 0;
}//END replaceIncumbent()

