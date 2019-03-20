/*
 * stocUpdate.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "stoc.h"

int stochasticUpdates(probType *prob, LPptr lp, basisType *basis, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, BOOL newOmegaFlag, int currentIter, double TOLERANCE, BOOL *newBasisFlag, BOOL subFeasFlag) {
	oneBasis *B;
	sparseVector dOmega;
	int 	cnt, lambdaIdx;
	BOOL	newSigmaFlag, newLambdaFlag, retainBasis;

	dOmega.cnt = prob->num->rvdOmCnt; dOmega.col = prob->coord->rvdOmCols;

	/* Update the column of delta structure if a new observation was encountered, and check the feasibility of existing bases with respect to new observation. */
	if ( newOmegaFlag ) {
		calcDelta(prob->num, prob->coord, lambda, delta, deltaRowLength, omega, newOmegaFlag, omegaIdx);

		/* Establish feasibility of basis with respect to current observations */
		dOmega.val = prob->coord->rvOffset[2]+omega->vals[omegaIdx];
		for ( cnt = 0; cnt < basis->cnt; cnt++ )
			basis->obsFeasible[cnt][omegaIdx] = checkBasisFeasibility(basis->vals[cnt], dOmega, prob->sp->senx, prob->num->cols, prob->num->rows, TOLERANCE);
	}

	if ( (B = newBasis(lp, prob->num->cols, prob->num->rows, currentIter, subFeasFlag)) == NULL ) {
		errMsg("algorithm", "stochasticUpdates", "failed to create a new basis type structure", 0);
		return -1;
	}

	/* check to see if the new basis was encountered before */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( B->feasFlag ) {
			if ( equalLongIntvec(B->cCode, basis->vals[cnt]->cCode, basis->cCodeLen) && equalLongIntvec(B->rCode,
					basis->vals[cnt]->rCode, basis->rCodeLen) ) {
				/* The basis is the same as one encountered before */
				freeOneBasis(B);
				basis->vals[cnt]->weight++;
				(*newBasisFlag) = FALSE;
#if defined (STOCH_CHECK)
				printf("An old basis encountered :: %d\n", cnt);
#endif
				return cnt;
			}
		}
	}

	if ( B->feasFlag ) {
		/* New basis encountered, fill the remainder of basis elements. */
		if ( prob->num->rvdOmCnt > 0 )
			calcBasis(lp, prob->num, prob->coord, prob->dBar, B, basis->basisDim);

		/* Decompose the dual solution into deterministic and stochastic components. */
		if ( decomposeDualSolution(lp, B, omega->vals[omegaIdx]+prob->coord->rvOffset[2], prob->num->rows) ) {
			errMsg("algorithm", "stochasticUpdates", "failed to decompose the dual solution", 0);
			return -1;
		}
	}
	else {
		if ( !(B->piDet = (vector) arr_alloc(prob->num->rows+1, double)) )
			errMsg("allocation", "decomposeDualSolution", "piS", 0);

		/* Record the dual and reduced cost on bounds. */
		if ( getDual(lp, B->piDet, prob->num->rows) ) {
			errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
			return 1;
		}
	}

	/* Elements of deterministic component of dual solution corresponding to rows with random elements in them */
	lambdaIdx = calcLambda(prob->num, prob->coord, B->piDet, lambda, &newLambdaFlag, TOLERANCE);

	/* Elements of deterministic component of dual solution with deterministic (mean value) right-hand side and transfer matrix. */
	B->sigmaIdx[0] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, B->piDet, B->mubBar,
			lambdaIdx, newLambdaFlag, currentIter, sigma, &newSigmaFlag, TOLERANCE);

	if ( newLambdaFlag )
		calcDelta(prob->num, prob->coord, lambda, delta, deltaRowLength, omega, FALSE, lambdaIdx);

	retainBasis = newSigmaFlag;
	for (cnt = 0; cnt < B->phiLength; cnt++ ) {
		/* Elements of basis column corresponding to rows with random elements in them */
		lambdaIdx = calcLambda(prob->num, prob->coord, B->phi[cnt], lambda, &newLambdaFlag, TOLERANCE);

		/* Compute the product of basis column with deterministic (mean value) right-hand side and transfer matrix. */
		B->sigmaIdx[cnt+1] = calcSigma(prob->num, prob->coord, prob->bBar, prob->Cbar, B->phi[cnt], 0,
				lambdaIdx, newLambdaFlag, currentIter, sigma, &newSigmaFlag, TOLERANCE);

		if ( newLambdaFlag )
			calcDelta(prob->num, prob->coord, lambda, delta, deltaRowLength, omega, FALSE, lambdaIdx);
		retainBasis = (retainBasis || newSigmaFlag);
	}

	if ( !retainBasis ) {
		/* All the sigmas computed were encountered before */
		for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
			if ( B->phiLength == basis->vals[cnt]->phiLength && basis->obsFeasible[cnt][omegaIdx] ) {
				if ( equalIntvec(B->sigmaIdx-1, basis->vals[cnt]->sigmaIdx-1, B->phiLength+1) ) {
					/* The basis was encountered before */
					freeOneBasis(B);
					basis->vals[cnt]->weight++;
					(*newBasisFlag) = FALSE;
					return cnt;
				}
			}
		}
	}

	/* Add the basis to the structure */
	basis->vals[basis->cnt] = B;

	if ( B->feasFlag ) {
		/* Establish feasibility of basis with respect to current observations */
		if ( !(basis->obsFeasible[basis->cnt] = (BOOL*) arr_alloc(deltaRowLength, BOOL)) )
			errMsg("allocation", "stochasticUpdates", "basis->obsFeasibility[n]", 0);
		for ( cnt = 0; cnt < omega->cnt; cnt++ ) {
			dOmega.val = prob->coord->rvOffset[2]+omega->vals[cnt];
			basis->obsFeasible[basis->cnt][cnt] = checkBasisFeasibility(B, dOmega, prob->sp->senx, prob->num->cols, prob->num->rows, TOLERANCE);
		}
	}
	else
		basis->obsFeasible[basis->cnt] = NULL;

	return basis->cnt++;

}//End stochasticUpdates()

/*This function loops through all the dual vectors found so far and returns the index of the one which satisfies the expression:
 * 				argmax { Pi x (R - T x X) | all Pi }
 * where X, R, and T are given.  It is calculated in this form:
 * 				Pi x bBar + Pi x bomega + (Pi x Cbar) x X + (Pi x Comega) x X.
 * Since the Pi's are stored in two different structures (sigma and delta), the index to the maximizing Pi is actually a structure
 * containing two indices.  (While both indices point to pieces of the dual vectors, sigma and delta may not be in sync with one
 * another due to elimination of non-distinct or redundant vectors. */
int computeIstar(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, vector piCbarX, vector Xvect, vector observ,
		int obs, int numSamples, BOOL pi_eval, double *argmax, BOOL isNew) {
	double 	arg, multiplier = 1.0;
	int 	cnt, maxCnt, c, basisUp, basisLow, sigmaIdx, lambdaIdx;

	if (pi_eval == TRUE)
		numSamples -= (int) (0.1*numSamples + 1);

	/* Establish the range of iterations over which the istar calculations are conducted. Only bases discovered in this iteration range are used. */
	if ( !isNew ) {
		basisUp = numSamples; basisLow = -INT_MAX;
	}
	else {
		basisUp = INT_MAX; basisLow = numSamples;
	}

	*argmax = -DBL_MAX; maxCnt = 0;

	/* Run through the list of basis to choose the one which provides the best lower bound */
	for ( cnt = 0; cnt < basis->cnt; cnt++ ) {
		if ( basis->vals[cnt]->feasFlag && basis->vals[cnt]->ck > basisLow && basis->vals[cnt]->ck <= basisUp ) {
			if ( basis->obsFeasible[cnt][obs] ) {
				arg = 0.0;
				for ( c = 0; c <= basis->vals[cnt]->phiLength; c++ ) {
					sigmaIdx = basis->vals[cnt]->sigmaIdx[c];
					lambdaIdx = sigma->lambdaIdx[sigmaIdx];
					if ( c == 0 )
						multiplier = 1.0;
					else
						multiplier = observ[coord->rvOffset[2] + basis->vals[cnt]->omegaIdx[c]];

					/* Start with (Pi x bBar) + (Pi x bomega) + (Pi x Cbar) x X */
					arg += multiplier*(sigma->vals[sigmaIdx].pib + delta->vals[lambdaIdx][obs].pib - piCbarX[sigmaIdx]);
					arg -= multiplier*vXv(delta->vals[lambdaIdx][obs].piC, Xvect, coord->rvCOmCols, num->rvCOmCnt);
				}

				if (arg > (*argmax)) {
					*argmax = arg;
					maxCnt = cnt;
				}
			}
		}
	}

	if ( (*argmax == -DBL_MAX ) )
		return -1;
	else
		return maxCnt;
}//END computeIstar

/* This function calculates a new column in the delta structure, based on a new observation of omega. Thus, lambda_pi X C and lambda_pi X b
 * are calculated for all values of lambda_pi, for the new C(omega) and b(omega).  Room in the array has already been allocated, so the function
 * only fills it, in the column specified by _obs_. It is assumed that this observation is distinct from all previous ones, and thus a new column
 * must be calculated. */
int calcDelta(numType *num, coordType *coord, lambdaType *lambda, deltaType *delta, int deltaRowLength, omegaType *omega, BOOL newOmegaFlag, int elemIdx) {
	sparseMatrix COmega;
	sparseVector bOmega;
	vector 		 lambdaPi, piCrossC;
	int 		 idx;

	/* extract the coordinates and number of random elements */
	bOmega.cnt = num->rvbOmCnt;	bOmega.col = coord->rvbOmRows;
	COmega.cnt = num->rvCOmCnt; COmega.col = coord->rvCOmCols; COmega.row = coord->rvCOmRows;

	if ( newOmegaFlag ) {
		/* Case I: New observation encountered. */
		bOmega.val= omega->vals[elemIdx];
		COmega.val = omega->vals[elemIdx] + num->rvbOmCnt;

		/* For all dual vectors, lambda(pi), calculate pi X bomega and pi X Comega */
		for (idx = 0; idx < lambda->cnt; idx++) {
			/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
			lambdaPi = expandVector(lambda->vals[idx], coord->rvRows, num->rvRowCnt, num->rows);

			/* Multiply the dual vector by the observation of bomega and Comega */
			/* Reduce PIxb from its full vector form into a sparse vector */
			delta->vals[idx][elemIdx].pib = vXvSparse(lambdaPi, &bOmega);
			if ( num->rvCOmCnt != 0 ) {
				piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
				delta->vals[idx][elemIdx].piC = reduceVector(piCrossC, coord->rvCOmCols, num->rvCOmCnt);
				mem_free(piCrossC);
			}
			else
				delta->vals[idx][elemIdx].piC = NULL;

			mem_free(lambdaPi);
		}
	}
	else {
		/* Case II: New dual vector encountered. */
		if ( !(delta->vals[elemIdx] = (pixbCType *) arr_alloc(deltaRowLength, pixbCType)))
			errMsg("allocation", "calcDeltaRow", "delta->val[cnt]", 0);

		/* expand the compressed lambda vector */
		lambdaPi = expandVector(lambda->vals[elemIdx], coord->rvRows, num->rvRowCnt, num->rows);

		/* go through all the observations and compute pi x b and pi x C */
		for (idx = 0; idx < omega->cnt; idx++) {

			bOmega.val= omega->vals[idx];
			COmega.val = omega->vals[idx] + num->rvbOmCnt;

			delta->vals[elemIdx][idx].pib = vXvSparse(lambdaPi, &bOmega);
			if ( num->rvCOmCnt != 0 ) {
				piCrossC = vxMSparse(lambdaPi, &COmega, num->prevCols);
				delta->vals[elemIdx][idx].piC = reduceVector(piCrossC, coord->rvCOmCols, num->rvCOmCnt);
				mem_free(piCrossC);
			}
			else
				delta->vals[elemIdx][idx].piC = NULL;
		}
		mem_free(lambdaPi);
	}

	return 0;
}//END calcDelta()

/* This function stores a new lambda_pi vector in the lambda structure.  Each lambda_pi represents only those dual variables whose rows in the
 * constraint matrix have random elements.  Thus  the (full) dual vector, Pi,  passed to the function is converted into the sparse vector lambda_pi.
 * This vector is then compared with all previous lambda_pi vectors, searching for a duplication. If a duplicate is found, the vector is not added
 * to the structure, and the function returns the index of the duplicate vector. Otherwise, it adds the vector to the end of the structure,
 *and returns an index to the last element in lambda. */
int calcLambda(numType *num, coordType *coord, vector Pi, lambdaType *lambda, BOOL *newLambdaFlag, double TOLERANCE) {
	int 	pi_idx;
	vector	lambda_pi;

	/* Pull out only those elements in dual vector which have rv's */
	lambda_pi = reduceVector(Pi, coord->rvRows, num->rvRowCnt);

	/* Compare resulting lambda_pi with all previous vectors */
	for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
		if (equalVector(lambda_pi, lambda->vals[pi_idx], num->rvRowCnt, TOLERANCE)) {
			mem_free(lambda_pi);
			*newLambdaFlag = FALSE;
			return pi_idx;
		}

	/* Add the vector to lambda structure */
	lambda->vals[lambda->cnt] = lambda_pi;
	*newLambdaFlag = TRUE;

	return lambda->cnt++;
}//END calcLambda

int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vector pi, double mubBar,
		int idxLambda, BOOL newLambdaFlag, int currentIter, sigmaType *sigma, BOOL *newSigmaFlag, double TOLERANCE) {
	vector	piCBar, temp;
	double 	pibBar;
	int 	cnt;

	/* sigma = \pi_t^\top \bar{b}_t - \bar{C}_t^\top \pi_t */
	pibBar = vXvSparse(pi, bBar) + mubBar;

	temp = vxMSparse(pi, CBar, num->prevCols);
	piCBar = reduceVector(temp, coord->CCols, num->cntCcols);
	mem_free(temp);

	if (!newLambdaFlag){
		for (cnt = 0; cnt < sigma->cnt; cnt++) {
			if (DBL_ABS(pibBar - sigma->vals[cnt].pib) <= TOLERANCE) {
				if (equalVector(piCBar, sigma->vals[cnt].piC, num->cntCcols, TOLERANCE))
					if(sigma->lambdaIdx[cnt]== idxLambda){
						mem_free(piCBar);
						(*newSigmaFlag) = FALSE;
						return cnt;
					}
			}
		}
	}

	(*newSigmaFlag) = TRUE;
	sigma->vals[sigma->cnt].pib  = pibBar;
	sigma->vals[sigma->cnt].piC  = piCBar;
	sigma->lambdaIdx[sigma->cnt] = idxLambda;
	sigma->ck[sigma->cnt] = currentIter;

	return sigma->cnt++;

}//END calcSigma()

/* This function obtains a new vector of realizations of the random variables. It compares the new vector with all previous vectors, looking for
 * a duplication.  If it finds a duplicate, it returns the index of that duplicate; otherwise, it adds the vector to the list of distinct realizations
 * and returns the index of that realization. Note that the simulated observation does not have contain one-norm, while the values stored in
 * omegaType do */
int calcOmega(vector observ, int begin, int end, omegaType *omega, BOOL *newOmegaFlag, double TOLERANCE) {
	int cnt;

	/* Compare vector with all the previous observations */
	for (cnt = 0; cnt < omega->cnt; cnt++)
		if (equalVector(observ, omega->vals[cnt], end-begin, TOLERANCE)) {
			(*newOmegaFlag) = FALSE;
			omega->weights[cnt]++;
			return cnt;
		}

	/* Add the realization vector to the list */
	omega->vals[omega->cnt] = duplicVector(observ, end-begin);
	omega->weights[omega->cnt] = 1;
	(*newOmegaFlag) = TRUE;

#ifdef STOCH_CHECK
	printf("Observation (%d): ", *newOmegaFlag);
	printVector(omega->vals[omega->cnt], end - begin, NULL);
#endif

	return omega->cnt++;
}//calcOmega()

/* This function compute the reduced cost of every second stage variables. They will be used to calculate the \mu x b and then added to the \pi x b. */
int computeMU(LPptr lp, intvec cstat, int numCols, double *mubBar) {
	vector	dj, u;
	int		n;

	(*mubBar) = 0.0;

	if ( !(dj = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "dual slacks", 0);
	if ( !(u = (vector) arr_alloc(numCols+1, double)))
		errMsg("allocation", "computeMu", "TDA solutions", 0);

	if ( getPrimal(lp, u, numCols) ) {
		errMsg("solver", "forOptPass", "failed to obtain primal solution", 0);
		return 1;
	}
	if (getDualSlacks(lp, dj, numCols) ) {
		errMsg("solver", "computeMu", "failed to obtain dual slacks", 0);
		return 1;
	}

	for (n = 1; n <= numCols;  n++) {
		switch (cstat[n]) {
		case AT_LOWER:
			(*mubBar) += dj[n]*u[n];
			break;
		case AT_UPPER:
			(*mubBar) += dj[n]*u[n];
			break;
		default:
			break;
		}
	}

	mem_free(u); mem_free(dj);

	return 0;
}//END compute_mu()

/* This function allocates a new lambda structure, with room for num_lambdas lambda vectors of size vect_size.  It returns a pointer to the structure.
 * Only some of the individual lambda vectors are expected to be allocated (according to the num_vect parameter) so that there is room for new
 * lambdas to be created. */
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows) {
	lambdaType *lambda;
	int cnt;

	if (!(lambda = (lambdaType *) mem_malloc (sizeof(lambdaType))))
		errMsg("allocation", "newLambda", "lambda",0);

	if (!(lambda->vals = arr_alloc(num_iter, vector)))
		errMsg("allocation", "newLambda", "lambda->val",0);

	for (cnt = 0; cnt < numLambda; cnt++)
		if (!(lambda->vals[cnt] = arr_alloc(numRVrows + 1, double)))
			errMsg("allocation", "newLambda", "lambda->val[cnt]",0);

	lambda->cnt = numLambda;

	return lambda;
}//END new_lambda

/* This function creates a new sigma structure, and allocates memory for the arrays associated with it.  It returns a pointer to this structure.
 * Some pi X T vectors are also allocated, according to the num_vals parameter  (num_vals is expected to be less than num_sigmas, so that there
 * is room for further work).  Note that  memory for sigma->col is not allocated, but is taken from prob.*/
sigmaType *newSigma(int numIter, int numNzCols, int numPi) {
	sigmaType *sigma;
	int cnt;

	if (!(sigma = (sigmaType *) mem_malloc (sizeof(sigmaType))))
		errMsg("allocation", "newSigma", "sigma",0);
	if (!(sigma->lambdaIdx = (intvec) arr_alloc(numIter, int)))
		errMsg("allocation", "newSigma", "sigma->lambIdx",0);
	if (!(sigma->ck = (intvec) arr_alloc(numIter, int)))
		errMsg("allocation", "newSigma", "sigma->ck",0);
	if (!(sigma->vals = arr_alloc(numIter, pixbCType)))
		errMsg("allocation", "newSigma", "sigma->vals",0);
	for (cnt = 0; cnt < numPi && cnt < numIter; cnt++)
		if (!(sigma->vals[cnt].piC = arr_alloc(numNzCols+1, double)))
			errMsg("allocation", "newSigma", "sigma->val[cnt]",0);

	sigma->cnt = numPi;

	return sigma;
}//END newSigma

/***********************************************************************\
 ** This function creates a new delta structure with arrays of the specified
 ** size and returns a pointer to it.  Note that the pi X T vectors
 ** themselves are not allocated, since they will not all be filled with
 ** values.  (they are only filled as they are produced).
 ** Not even the arrays of pi_R_T_types are allocated, as this also
 ** occurs in calc_delta_row().  However, the column coordinates of the
 ** eventual multiplications are initialized, since they are known.
 \***********************************************************************/
deltaType *newDelta(int numIter) {
	deltaType *delta;

	if (!(delta = (deltaType *) mem_malloc (sizeof(deltaType))))
		errMsg("Allocation", "newDelta", "d",0);
	if (!(delta->vals = (pixbCType **) arr_alloc(numIter, pixbCType *)))
		errMsg("Allocation", "newDelta", "d->val",0);

	return delta;
}//END newDelta

/* This function allocates memory for an omega structure.  It allocates the memory to structure elements: a vector to hold an array of
 * observation and the weights associated with it. */
omegaType *newOmega(int numOmega, int numIter) {
	omegaType *omega;

	if ( !(omega = (omegaType *) mem_malloc(sizeof(omegaType))) )
		errMsg("allocation","newOmega", "omega", 0);
	if ( !(omega->weights = (intvec) arr_alloc(numIter, int)) )
		errMsg("allocation", "newOmega", "omega->weights", 0);
	if ( !(omega->vals = (vector *) arr_alloc(numIter, vector)) )
		errMsg("allocation", "newOmega", "omega->vals", 0);
	omega->numRV = numOmega;
	omega->cnt = 0;

	return omega;
}//END newOmega()

void freeOmegaType(omegaType *omega, BOOL partial) {
	int n;

	if ( omega->vals ) {
		for ( n = 0; n < omega->cnt; n++ )
			if ( omega->vals[n] ) mem_free(omega->vals[n]);
		if ( partial ) {
			omega->cnt = 0;
			return;
		}
		mem_free(omega->vals);
	}
	if ( omega->weights ) mem_free(omega->weights);
	mem_free(omega);

}//END freeOmegaType()

void freeLambdaType(lambdaType *lambda, BOOL partial) {
	int n;

	if (lambda) {
		if (lambda->vals) {
			for ( n = 0; n < lambda->cnt; n++ )
				if (lambda->vals[n]) mem_free(lambda->vals[n]);
			if ( partial ) {
				lambda->cnt = 0;
				return;
			}
			mem_free(lambda->vals);
		}
		mem_free(lambda);
	}

}//END freeLambdaType()

void freeSigmaType(sigmaType *sigma, BOOL partial) {
	int n;

	if (sigma) {
		for ( n = 0; n < sigma->cnt; n++ )
			if (sigma->vals[n].piC) mem_free(sigma->vals[n].piC);
		if ( partial ) {
			sigma->cnt = 0;
			return;
		}
		if (sigma->lambdaIdx) mem_free(sigma->lambdaIdx);
		if (sigma->vals) mem_free(sigma->vals);
		if (sigma->ck) mem_free(sigma->ck);
		mem_free(sigma);
	}

}//END freeSigmaType()

void freeDeltaType (deltaType *delta, int numDeltaRows, int omegaCnt, BOOL partial) {
	int n, m;

	if (delta) {
		if (delta->vals) {
			for ( n = 0; n < numDeltaRows; n++ ) {
				if (delta->vals[n]) {
					for ( m = 0; m < omegaCnt; m++ )
						if (delta->vals[n][m].piC)
							mem_free(delta->vals[n][m].piC);
					mem_free(delta->vals[n]);
				}
			}
			if ( partial )
				return;
			mem_free(delta->vals);
		}
		mem_free(delta);
	}

}//END freeDeltaType()
