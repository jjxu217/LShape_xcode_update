/*
 * master.c
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

/* This function is the regularized QP version of master problem. The master problem is solved after the newest cut is added to master problem,
 the incumbent cut is updated if necessary. Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, double lb) {
	double 	d2 = 0.0; /* height at the candidate solution. */
	int 	status, i;

	if( changeEtaCol(cell->master->lp, num->rows, num->cols, cell->k, cell->cuts) ) {
		errMsg("algorithm", "solveQPMaster", "failed to change the eta column coefficients", 0);
		return 1;
	}

	if ( cell->lbType == NONTRIVIAL ) {
		/* update the right-hand side of cuts to reflect the non-trivial lower bound */
		if ( updateRHS(cell->master->lp, cell->cuts, cell->k, cell->lb) ) {
			errMsg("algorithm", "solveQPMaster", "failed to update right-hand side with lower bound information", 0);
			return 1;
		}
	}

#ifdef ALGO_CHECK
	writeProblem(cell->master->lp, "cellMaster.lp");
#endif

	/* solve the master problem */
	clock_t tic = clock();
	changeQPSolverType(ALG_CONCURRENT);
	if ( solveProblem(cell->master->lp, cell->master->name, config.MASTER_TYPE, &status) ) {
		if ( status == STAT_INFEASIBLE ) {
			errMsg("algorithm", "solveQPMaster", "Master problem is infeasible. Check the problem formulation!",0);
			writeProblem(cell->master->lp, "infeasibleM.lp");
		}
		else {
			writeProblem(cell->master->lp, "errorM.lp");
			errMsg("algorithm", "solveQPMaster", "failed to solve the master problem", 0);
		}
		return 1;
	}
	cell->time.masterIter = ((double) (clock() - tic))/CLOCKS_PER_SEC;

	/* Get the most recent optimal solution to master program */
	if ( getPrimal(cell->master->lp, cell->candidX, num->cols) ) {
		errMsg("algorithm", "solveQPMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* add the incumbent back to change from \Delta X to X */
	for (i = 1; i <= num->cols; i++)
		d2 += cell->candidX[i] * cell->candidX[i];
	addVectors(cell->candidX, cell->incumbX, NULL, num->cols);

	/* update d_norm_k in soln_type. */
	if (cell->k == 1)
		cell->normDk_1 = d2;
	cell->normDk = d2;

	/* Get the dual solution too */
	if ( getDual(cell->master->lp, cell->piM, cell->master->mar) ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
		return 1;
	}

	if ( getDualSlacks(cell->master->lp, cell->djM, num->cols) ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta(xbar + \Delta X) */
	cell->candidEst = vXvSparse(cell->candidX, dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, num->cols, lb);

	/* Calculate gamma for next improvement check on incumbent x. */
	cell->gamma = cell->candidEst - cell->incumbEst;

	return 0;
}//END solveQPMaster()

int addCut2Master(oneProblem *master, oneCut *cut, vector vectX, int lenX) {
	intvec 	indices;
	int 	cnt;
	static int cummCutNum = 0;

	/* Set up indices */
	if (!(indices = (intvec) arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

	/* Cut right-hand side */
	if ( config.MASTER_TYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, vectX, NULL, lenX);

	/* Set up the cut name */
	sprintf(cut->name, "cut_%04d", cummCutNum++);

	/* add the cut to the cell cuts structure as well as on the solver */
	if ( addRow(master->lp, lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta, cut->name) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return 1;
	}
	cut->rowNum = master->mar++;

#ifdef CUT_CHECK
	writeProblem(cell->master->lp,"master_wNewCut.lp");
#endif

	mem_free(indices);
	return 0;
}//END addCuts2Master()

int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar) {

	if ( changeQPproximal(cell->master->lp, prob->num->cols, quadScalar) ) {
		errMsg("algorithm", "algoIntSD", "failed to change the proximal term", 0);
		return 1;
	}

	if ( changeQPrhs(prob, cell, incumbX) ) {
		errMsg("algorithm", "algoIntSD", "failed to change the right-hand side to convert the problem into QP", 0);
		return 1;
	}

	if ( changeQPbds(cell->master->lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX, 0) ) {
		errMsg("algorithm", "algoIntSD", "failed to change the bounds to convert the problem into QP", 0);
		return 1;
	}

	return 0;
}//END constructQP()

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts) {
	double	coef[1];
	int 	c;

	for (c = 0; c < cuts->cnt; c++){
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[0] = (double) (k) / (double) cuts->vals[c]->numSamples;         // coefficient k/j of eta column

		if ( changeCol(lp, numCols, coef, cuts->vals[c]->rowNum, cuts->vals[c]->rowNum+1) ) {
			errMsg("solver", "changeEtaCol", "failed to change eta column in the stage problem", 0);
			return 1;
		}
	}

	return 0;
}//END changeEtaCol()

int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb) {
	int 	cnt;
	vector	rhs;
	intvec	indices;

	if (!(rhs = arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);
	if (!(indices = arr_alloc(cuts->cnt, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		rhs[cnt] = cuts->vals[cnt]->alphaIncumb + ((double) numIter / (double) cuts->vals[cnt]->numSamples - 1) * lb;
		indices[cnt] = cuts->vals[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	if ( changeRHS(lp, cuts->cnt, indices, rhs) ) {
		errMsg("solver", "updateRHS", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int changeQPproximal(LPptr lp, int numCols, double sigma) {
	int    n;
	vector qsepvec;

	if (!(qsepvec = arr_alloc(numCols+1, double)))
		errMsg("Allocation", "changeQPproximal", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (n = 0; n < numCols; n++)
		qsepvec[n] = 0.5 * sigma;
	qsepvec[n] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	if ( copyQPseparable(lp, qsepvec) ) {
		errMsg("solver", "changeQPproximal", "failed to copy Q matrix", 0);
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END constructQPproximal

/* In the regularized QP method, we need to change the rhs of x to d. The
 * 		 A * x 			= b
 * 		 eta + beta * x >= alpha
 * Since x = xbar + d, the corresponding changes will therefore be:
 * 		 A * d = b - A * xbar
 * 		 eta + beta * d >= alpha - beta * xbar
 * But as long as the incumbent sulotion does not change, b - A * xbar and alpha - beta * xbar (for the existing cuts) won't change. So we only need
 * to change it when the incumbent changes.
 *
 * On the other hand, in each iteration, a new cut will be added (and/or some cuts may be dropped) and therefore we need to shift the rhs of the
 * added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of in the routine addCut() in cuts.c. We do not need to worry about the shift
 * of rhs for the dropped cuts.
 * This function performs the change of rhs when the incumbent changes, as described above. */
int changeQPrhs(probType *prob, cellType *cell, vector xk) {
	int 	status = 0, cnt;
	vector 	rhs;
	intvec 	indices;

	if (!(rhs =(vector) arr_alloc(prob->num->rows+cell->cuts->cnt+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(intvec) arr_alloc(prob->num->rows+cell->cuts->cnt, int)))
		errMsg("Allocation", "changeRhs", "indices",0);
	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
	 vector, and the x vector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, xk, rhs);

	/*** new rhs = alpha - beta * xbar (benders cuts)***/
	for (cnt = 0; cnt < cell->cuts->cnt; cnt++) {
		rhs[prob->num->rows+cnt+1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, xk, NULL, prob->sp->mac);
		indices[prob->num->rows+cnt] = cell->cuts->vals[cnt]->rowNum;

		cell->cuts->vals[cnt]->alphaIncumb = rhs[prob->num->rows+cnt+1];
	}

	/* Now we change the right-hand of the master problem. */
	status = changeRHS(cell->master->lp, prob->num->rows + cell->cuts->cnt, indices, rhs+1);
	if (status)	{
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrhs()

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk, int offset) {
	int 	status = 0, cnt;
	vector	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt+offset;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt+offset;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
oneProblem *newMaster(oneProblem *orig, double lb) {
	oneProblem 	*master;
	int         r, i, j, idx, cnt;
	long        colOffset, rowOffset;
	char        *q;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to master -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	master->type 	= config.MASTER_TYPE;               /* type of problem: LP, QP, MIP or MIQP */
	master->objsen 	= orig->objsen;                 	/* sense of the objective: 1 for minimization and -1 for maximization */
	master->mar 	= orig->mar;                       	/* number of rows */
	master->numInt 	= orig->numInt;                 	/* number of integer variables in the problem  */
	master->numnz 	= orig->numnz;                   	/* number of non-zero elements in constraint matrix */
	master->matsz 	= orig->matsz;                   	/* extended matrix size */
	master->marsz 	= orig->marsz;                   	/* extended row size */
	master->rstorsz = orig->rstorsz;               		/* memory size for storing row names */
	master->mac 	= orig->mac+1;           			/* number of columns + etas */
	master->macsz 	= orig->macsz + 1;       			/* extended column size */
	master->cstorsz 	= orig->cstorsz + NAMESIZE;    	/* memory size for storing column names */

	/* Allocate memory to the information whose type is string */
	if (!(master->name = (string) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->name",0);
	if (!(master->senx = (string) arr_alloc(master->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->senx",0);
	if (!(master->ctype = (string) arr_alloc(master->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->ctype",0);
	if (!(master->objname = (string) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objname",0);
	if (!(master->cname = (string*) arr_alloc(master->macsz,string)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cname",0);
	if (!(master->cstore = (string) arr_alloc(master->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cstore",0);
	if ( master->mar > 0 ) {
		if (!(master->rname = (string *) arr_alloc(master->marsz,string)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rname",0);
		if (!(master->rstore = (string) arr_alloc(master->rstorsz, char)))
			errMsg("Allocation", "new_master", "Fail to allocate memory to master->rstore",0);
	}
	else {
		master->rname = NULL; master->rstore = NULL;
	}

	/* Allocate memory to the information whose type is vector */
	if (!(master->objx = (vector) arr_alloc(master->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objx",0);
	if (!(master->rhsx = (vector) arr_alloc(master->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rhsx",0);
	if (!(master->matval = (vector) arr_alloc(master->matsz, double)))
		errMsg("allocation", "new_master", "master->matval",0);
	if (!(master->bdl = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdl",0);
	if (!(master->bdu = (vector) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdu",0);

	/* Allocate memory to the information whose type is intvec */
	if (!(master->matbeg = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matbeg",0);
	if (!(master->matcnt = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matcnt",0);
	if (!(master->matind = (intvec) arr_alloc(master->matsz, int)))
		errMsg("allocation", "new_master", "master->matind",0);

	strcpy(master->name, orig->name);           /* Copy problem name */
	strcpy(master->objname, orig->objname);     /* Copy objective name */

	/* Copy problem's column and row names. Calculate difference in pointers for master/copy row and column names. */
	i = 0;
	for (q = orig->cname[0]; q < orig->cname[0] + orig->cstorsz; q++)
		master->cstore[i++] = *q;
	colOffset = master->cstore - orig->cname[0];

	if ( master->mar > 0 ) {
		i = 0;
		for (q = orig->rname[0]; q < orig->rname[0] + orig->rstorsz; q++)
			master->rstore[i++] = *q;
		rowOffset = master->rstore - orig->rname[0];
	}

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < orig->mac; j++) {
		master->objx[j] = orig->objx[j];		/* Copy objective function coefficients */
		master->ctype[j] = orig->ctype[j];		/* Copy the decision variable type */
		master->bdu[j] = orig->bdu[j];			/* Copy the upper bound and lower bound */
		master->bdl[j] = orig->bdl[j];
		master->cname[j] = orig->cname[j] + colOffset; /* Copy column names, offset by length */
		master->matbeg[j] = cnt;				/* Copy the master sparse matrix beginning position of each column */
		master->matcnt[j] = orig->matcnt[j];	/* Copy the sparse matrix non-zero element count */
		master->ctype[j] = orig->ctype[j];
		/* Loop through all non-zero elements in this column */
		for (idx = orig->matbeg[j]; idx < orig->matbeg[j] + orig->matcnt[j]; idx++) {
			master->matval[cnt] = orig->matval[idx];	/* Copy the non-zero coefficient */
			master->matind[cnt] = orig->matind[idx];	/* Copy the row entry of the non-zero elements */
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < orig->mar; r++) {
		master->rhsx[r] = orig->rhsx[r];		/* Copy the right hand side value */
		master->senx[r] = orig->senx[r];		/* Copy the constraint sense */
		master->rname[r] = orig->rname[r] + rowOffset;	/* Copy row names, offset by length */
	}

	/* Initialize information for the extra column in the new master. */
	colOffset = orig->cstorsz;
	strcpy(master->cstore + orig->cstorsz, "eta");
	master->cname[orig->mac] = master->cstore + colOffset;
	master->objx[orig->mac] = 1.0;			// orig->mac is the last column in the original master
	master->ctype[orig->mac] = 'C';
	master->bdu[orig->mac] = INFBOUND;
	master->bdl[orig->mac] = lb;
	master->matbeg[orig->mac] = orig->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
	master->matcnt[orig->mac] = 0;               // Only optimality cuts has eta

	/* Load the copy into CPLEX */
	master->lp = setupProblem(master->name, master->type, master->mac, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if defined(SETUP_CHECK)
	if ( writeProblem(master->lp, "newMaster.lp") ) {
		errMsg("solver", "newMaster", "failed to write master problem to file", 0);
		return NULL;
	}
#endif

	return master;

}//END newMaster()
