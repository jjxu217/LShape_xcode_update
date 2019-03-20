/*
 * subprob.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "stoc.h"

/* This function will solve a new subproblem. This involves replacing the right-hand side of the subproblem with new values, based upon some
 * observation of omega, and some X vector of primal variables from the master problem.  Generally, the latest observation is used.  When
 * forming a normal cut, the candidate x should be used, while the incumbent x should be used for updating the incumbent cut. */
int solveSubprob(probType *prob, oneProblem *subproblem, vector Xvect, basisType *basis, lambdaType *lambda, sigmaType *sigma, deltaType *delta, int deltaRowLength,
		omegaType *omega, int omegaIdx, BOOL *newOmegaFlag, int currentIter, double TOLERANCE, BOOL *subFeasFlag, BOOL *newBasisFlag,
		double *subprobTime, double *argmaxTime) {
	int  	status;
	clock_t tic;

	/* (a) compute and change the right-hand side using current observation and first-stage solution */
	if ( computeRHS(subproblem->lp, prob->num, prob->coord, prob->bBar, prob->Cbar, Xvect, omega->vals[omegaIdx]) ) {
		errMsg("algorithm", "solveSubprob", "failed to compute subproblem right-hand side", 0);
		return 1;
	}

	if ( prob->num->rvdOmCnt > 0 ) {
		/* (b) Compute and change the cost coefficients using current observation */
		if ( computeCostCoeff(subproblem->lp, prob->num, prob->coord, prob->dBar, omega->vals[omegaIdx]) ) {
			errMsg("algorithm", "solveSubprob", "failed to compute subproblem cost coefficients", 0);
			return 1;
		}
	}

#if defined(ALGO_CHECK)
	writeProblem(subproblem->lp, "subproblem.lp");
#endif

	/* (c) Solve the subproblem to obtain the optimal dual solution. */
	tic = clock();
	setIntParam(PARAM_PREIND, OFF);
	changeLPSolverType(ALG_PRIMAL);
	if ( solveProblem(subproblem->lp, subproblem->name, subproblem->type, &status) ) {
		if ( status == STAT_INFEASIBLE ) {
			/* Set the subproblem feasibility flag to false and proceed to complete stochastic updates. These updates are
			 * used to generate the feasibility cuts later. */
			printf("Subproblem is infeasible for current first-stage decision and observation.\n");
			writeProblem(subproblem->lp, "infeasibleSP.lp");
			(*subFeasFlag) = FALSE;
		}
		else {
			errMsg("algorithm", "solveSubprob", "failed to solve subproblem in solver", 0);
			return 1;
		}
	}
	setIntParam(PARAM_PREIND, ON);
	(*subprobTime) += ((double) (clock() - tic))/CLOCKS_PER_SEC;

#ifdef STOCH_CHECK
	double obj;
	obj = getObjective(subproblem->lp, PROB_LP);
	printf("Objective value of Subproblem  = %lf\n", obj);
#endif

	if ( newBasisFlag!= NULL ) {
		tic = clock();
		/* (d) update the stochastic elements in the problem */
		status = stochasticUpdates(prob, subproblem->lp, basis, lambda, sigma, delta, deltaRowLength,
				omega, omegaIdx, (*newOmegaFlag), currentIter, TOLERANCE ,newBasisFlag, (*subFeasFlag));
		(*newOmegaFlag) = FALSE;
		(*argmaxTime) += ((double) (clock()-tic))/CLOCKS_PER_SEC;

#ifdef STOCH_CHECK
		obj = sigma->vals[status].pib - vXv(sigma->vals[status].piC, Xvect, prob->coord->CCols, prob->num->cntCcols);
		obj += delta->vals[sigma->lambdaIdx[status]][omegaIdx].pib - vXv(delta->vals[sigma->lambdaIdx[status]][omegaIdx].piC,
				omega->vals[omegaIdx], prob->coord->rvCOmCols, prob->num->rvCOmCnt);
		printf("Objective function estimate    = %lf\n", obj);
#endif
	}

	return 0;
}// END solveSubprob()

/* This function computes the right hand side of the subproblem, based on a given X vector and a given observation of omega.
 * It is defined as:
 * 			rhs = R(omega) - T(omega) x X
 * and is calculated as:
 * 			rhs = (Rbar - Tbar x X) + (Romega - Tomega x X)
 *
 * where the "bar" denotes the fixed or mean value, and the "omega" denotes a random variation from this mean. The function allocates an array
 * for the vector, which must be freed by the customer.  Also, the zeroth position of this rhs vector is reserved, and the actual values begin at rhs[1].
 * R is b, and T is C
 \***********************************************************************/
int computeRHS(LPptr lp, numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vector X, vector obs) {
	sparseMatrix Comega;
	sparseVector bomega;
	vector rhs;
	int cnt, *indices;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)) )
		errMsg("allocation", "solveSubprob", "indices", 0);
	for ( cnt = 0; cnt < num->rows; cnt++ )
		indices[cnt] = cnt;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows; bomega.val = obs + coord->rvOffset[0];

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols + num->rvbOmCnt;
	Comega.row = coord->rvCOmRows + num->rvbOmCnt; Comega.val = obs + coord->rvOffset[1];

	/* Start with the values of b(omega) -- both fixed and varying */
	rhs = expandVector(bBar->val, bBar->col, bBar->cnt, num->rows);
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* (cumulatively) subtract values of C(omega) x X -- both fixed and varying */
	rhs = MSparsexvSub(Cbar, X, rhs);
	rhs = MSparsexvSub(&Comega, X, rhs);

	/* change the right-hand side in the solver */
	if ( changeRHS(lp, num->rows, indices, rhs + 1) ) {
		errMsg("solver", "solveSubprob", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	mem_free(indices); mem_free(rhs);
	return 0;
}//END computeRHS()

int computeCostCoeff(LPptr lp, numType *num, coordType *coord, sparseVector *dBar, vector observ) {
	sparseVector dOmega;
	vector cost;
	int	cnt, *indices;

	if ( !(indices = (intvec) arr_alloc(num->cols, int)) )
		errMsg("allocation", "solveSubprob", "indices", 0);
	for ( cnt = 0; cnt < num->cols; cnt++ )
		indices[cnt] = cnt;

	dOmega.cnt = num->rvdOmCnt; dOmega.col = coord->rvdOmCols; dOmega.val = coord->rvOffset[2] + observ;

	/* Extract the cost coefficients */
	cost = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
	for (cnt = 1; cnt <= dOmega.cnt; cnt++)
		cost[dOmega.col[cnt]] += dOmega.val[cnt];

	/* change cost coefficients in the solver */
	if ( changeObjx(lp, num->cols, indices, cost+1) ) {
		errMsg("solver", "solve_subprob", "failed to change the cost coefficients in the solver",0);
		return -1;
	}

	mem_free(indices); mem_free(cost);
	return 0;
}//END computeCostCoeff()

void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vector rhs, vector X) {
	int cnt;

	/* copy the original right-hand side */
	for (cnt = 1; cnt <= bBar->cnt; cnt++)
		rhs[bBar->col[cnt]] = bBar->val[cnt];

	/* change the right-hand side with first stage solution */
	rhs = MSparsexvSub(Cbar, X, rhs);

}//END chgRHSwMean()

int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vector observ, vector spRHS, vector X) {
	sparseVector bomega;
	sparseMatrix Comega;
	vector 	rhs;
	intvec	indices;
	int		cnt, stat1;

	bomega.cnt = num->rvbOmCnt;	bomega.col = coord->rvbOmRows; bomega.val = observ;

	Comega.cnt = num->rvCOmCnt; Comega.col = coord->rvCOmCols + num->rvbOmCnt;
	Comega.row = coord->rvCOmRows + num->rvbOmCnt; Comega.val = observ + num->rvbOmCnt;

	if ( !(indices = (intvec) arr_alloc(num->rows, int)) )
		errMsg("allocation", "chgRHSwObserv", "indices", 0);
	if ( !(rhs = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "chgRHSwObserv", "rhs", 0);

	/* copy right-hand side modified with mean information */
	for ( cnt = 1; cnt <= num->rows; cnt++ ) {
		rhs[cnt] = spRHS[cnt];
		indices[cnt-1] = cnt-1;
	}

	/* change right-hand side with randomness in b */
	for (cnt = 1; cnt <= bomega.cnt; cnt++)
		rhs[bomega.col[cnt]] += bomega.val[cnt];

	/* change right-hand side with randomness in transfer matrix */
	rhs = MSparsexvSub(&Comega, X, rhs);

	/* change the right-hand side in the solver */
	stat1 = changeRHS(lp, num->rows, indices, rhs + 1);
	if ( stat1 ) {
		errMsg("solver", "chgRHSwObserv", "failed to change the right-hand side in the solver",0);
		return 1;
	}

	mem_free(rhs); mem_free(indices);
	return 0;

}//END chgRHSwRand()

int chgObjxwObserv(LPptr lp, numType *num, coordType *coord, vector cost, intvec indices, vector observ) {
	vector vals;
	int n;

	if ( !(vals = (vector) arr_alloc(num->rvdOmCnt+1, double)) )
		errMsg("allocation", "chgObjwObserv", "vals", 0);

	for ( n = 1; n <= num->rvdOmCnt; n++ )
		vals[n] = cost[n] + observ[coord->rvOffset[2]+n];

	if ( changeObjx(lp, num->rvdOmCnt, indices+1, vals+1) ) {
		errMsg("solver", "chgObjswObserv", "failed to change the cost coefficients in the solver",0);
		return 1;
	}

	mem_free(vals);
	return 0;
}//END chgObjwObserv()

oneProblem *newSubprob(oneProblem *sp) {

	/* since the basic structure of subproblem is not modified during the course of the algorithm, we just load it onto the solver */
	sp->lp = setupProblem(sp->name, sp->type, sp->mac, sp->mar, sp->objsen, sp->objx, sp->rhsx, sp->senx,sp->matbeg, sp->matcnt, sp->matind, sp->matval, sp->bdl, sp->bdu, NULL, sp->cname, sp->rname, sp->ctype);
	if ( sp->lp == NULL ) {
		errMsg("Problem Setup", "newSubprob", "sp",0);
		return NULL;
	}

#if 0
	int     status;
	char probName[NAMESIZE];
	sprintf(probName,"newSubprob%d.lp", agent);
	status = writeProblem(scell->sp->lp, probName);
	if ( status ) {
		errMsg("write problem", "new_subprob", "failed to write subproblems problem to file",0);
		return NULL;
	}
#endif

	return sp;
}//END new_subprob
