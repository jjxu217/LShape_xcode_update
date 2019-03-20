/*
 * algo.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern string outputDir;
extern configType config;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, string inputDir, string probName) {
	vector	 meanSol = NULL;
	probType **prob = NULL;
	cellType *cell = NULL;
	batchSummary *batch = NULL;
	FILE 	*sFile = NULL, *iFile = NULL;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol) )
		goto TERMINATE;

	printf("Starting two-stage stochastic decomposition.\n");
	sFile = openFile(outputDir, "results.dat", "w");
	iFile = openFile(outputDir, "incumb.dat", "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for ( int rep = 0; rep < config.NUM_REPS; rep++ ) {
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "Replication-%d\n", rep+1);
		fprintf(stdout, "\n====================================================================================================================================\n");
		fprintf(stdout, "Replication-%d\n", rep+1);

		/* setup the seed to be used in the current iteration */
		config.RUN_SEED[0] = config.RUN_SEED[rep+1];
		config.EVAL_SEED[0] = config.EVAL_SEED[rep+1];

		if ( rep != 0 ) {
			/* clean up the cell for the next replication */
			if ( cleanCellType(cell, prob[0], meanSol) ) {
				errMsg("algorithm", "algo", "failed clean the problem cell", 0);
				goto TERMINATE;
			}
		}

		clock_t tic = clock();
		/* Use two-stage stochastic decomposition algorithm to solve the problem */
		if ( solveCell(stoc, prob, cell) ) {
			errMsg("algorithm", "algo", "failed to solve the cell using 2-SD algorithm", 0);
			goto TERMINATE;
		}
		cell->time.repTime = ((double) clock() - tic)/CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		if (rep == 0 ) {
			writeOptimizationSummary(sFile, iFile, prob, cell, TRUE);
			writeOptimizationSummary(stdout, NULL, prob, cell, TRUE);
		}
		else {
			writeOptimizationSummary(sFile, iFile, prob, cell, FALSE);
			writeOptimizationSummary(stdout, NULL, prob, cell, FALSE);
		}

		/* evaluate the optimal solution*/
		if (config.EVAL_FLAG == 1)
			evaluate(sFile, stoc, prob, cell->subprob, cell->incumbX);

		/* Save the batch details and build the compromise problem. */
		if ( config.MULTIPLE_REP ) {
			buildCompromise(prob[0], cell, batch);
		}
	}

	if ( config.MULTIPLE_REP ) {
		/* Solve the compromise problem. */
		if ( solveCompromise(prob[0], batch)) {
			errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
			goto TERMINATE;
		}

		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		fprintf(sFile, "\n====================================================================================================================================\n");
		fprintf(sFile, "\n----------------------------------------- Compromise solution --------------------------------------\n\n");
		/* Evaluate the compromise solution */
		evaluate(sFile, stoc, prob, cell->subprob, batch->compromiseX);

		fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
		/* Evaluate the average solution */
		evaluate(sFile, stoc, prob, cell->subprob, batch->avgX);
	}

	fclose(sFile); fclose(iFile);
	printf("\nSuccessfully completed two-stage stochastic decomposition algorithm.\n");

	/* free up memory before leaving */
	if (meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

	TERMINATE:
	if(meanSol) mem_free(meanSol);
	freeBatchType(batch);
	freeCellType(cell);
	freeProbType(prob, 2);
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell) {
	vector 	observ;
	int		m, omegaIdx, candidCut;
	BOOL 	newOmegaFlag;
	clock_t	tic;

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "solveCell", "observ", 0);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == FALSE && cell->k < config.MAX_ITER) {
		cell->k++;
		tic = clock();
#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\nIteration-%d :: \n", cell->k);
#else
		if ( (cell->k -1) % 100 == 0) {
			printf("\nIteration-%4d: ", cell->k);
		}
#endif

		/******* 1. Optimality tests *******/
		if (optimal(prob, cell))
			break;

		/******* 2. Generate new observation, and add it to the set of observations *******/
		/* (a) Use the stoc file to generate observations */
		generateOmega(stoc, observ, config.TOLERANCE, &config.RUN_SEED[0]);

		/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
		for ( m = 0; m < stoc->numOmega; m++ )
			observ[m] -= stoc->mean[m];

		/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
		omegaIdx = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &newOmegaFlag, config.TOLERANCE);

		/******* 3. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formSDCut(prob, cell, cell->candidX, omegaIdx, &newOmegaFlag, prob[0]->lb)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			goto TERMINATE;
		}

		/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
		if (((cell->k - cell->iCutUpdt) % config.TAU == 0 ) ) {
			if ( (cell->iCutIdx = formSDCut(prob, cell, cell->incumbX, omegaIdx, &newOmegaFlag, prob[0]->lb) ) < 0 ) {
				errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
				goto TERMINATE;
			}
			cell->iCutUpdt = cell->k;
		}

		/******* 5. Check improvement in predicted values at candidate solution *******/
		if ( !(cell->incumbChg) && cell->k > 1)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement(prob[0], cell, candidCut);

		/******* 6. Solve the master problem to obtain the new candidate solution */
		if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->lb) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			goto TERMINATE;
		}

		cell->time.masterAccumTime += cell->time.masterIter; cell->time.subprobAccumTime += cell->time.subprobIter;
		cell->time.argmaxAccumTime += cell->time.argmaxIter; cell->time.optTestAccumTime += cell->time.optTestIter;
		cell->time.masterIter = cell->time.subprobIter = cell->time.optTestIter = cell->time.argmaxIter = 0.0;
		cell->time.iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time.iterAccumTime += cell->time.iterTime;
	}//END while loop

	mem_free(observ);
	return 0;

	TERMINATE:
	mem_free(observ);
	return 1;
}//END solveCell()

void writeOptimizationSummary(FILE *soln, FILE *incumb, probType **prob, cellType *cell, BOOL header) {

	if ( header ) {
		fprintf(soln, "\n--------------------------------------- Problem Information ----------------------------------------\n\n");
		fprintf(soln, "Problem                                : %s\n", prob[0]->name);
		fprintf(soln, "First Stage Rows                       : %d\n", prob[0]->num->rows);
		fprintf(soln, "First Stage Columns                    : %d\n", prob[0]->num->cols);
	}

	fprintf(soln, "\n------------------------------------------- Optimization -------------------------------------------\n\n");

	fprintf(soln, "Algorithm                              : Two-stage Stochastic Decomposition\n");
	fprintf(soln, "Number of iterations                   : %d\n", cell->k);
	fprintf(soln, "Lower bound estimate                   : %f\n", cell->incumbEst);
	fprintf(soln, "Total time                             : %f\n", cell->time.repTime);
	fprintf(soln, "Total time to solve master             : %f\n", cell->time.masterAccumTime);
	fprintf(soln, "Total time to solve subproblems        : %f\n", cell->time.subprobAccumTime);
	fprintf(soln, "Total time in argmax procedure         : %f\n", cell->time.argmaxAccumTime);
	fprintf(soln, "Total time in verifying optimality     : %f\n", cell->time.optTestAccumTime);

	if ( incumb != NULL ) {
		printVector(cell->incumbX, prob[0]->num->cols, incumb);
	}

}//END WriteStat



