/*
 * algo.c
 *
 *  Created on: Sep 21, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "benders.h"

extern string outputDir;
extern configType config;

#if defined(SAVE_DUALS)
dualsType *duals = NULL;
#endif

int algo (oneProblem *orig, timeType *tim, stocType *stoc, string probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	vector 	 meanSol;
	int 	 rep, m, n;
	FILE 	*sFile, *iFile = NULL;
	clock_t	tic;

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &meanSol) )
		goto TERMINATE;

	printf("Starting Benders decomposition.\n");
	sFile = openFile(outputDir, "results.dat", "w");
	if ( config.MASTER_TYPE == PROB_QP )
		iFile = openFile(outputDir, "incumb.dat", "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);

	for ( rep = 0; rep < config.NUM_REPS; rep++ ) {
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
				errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
				goto TERMINATE;
			}
		}

		/* Update omega structure */
		if ( config.SAA ) {
			setupSAA(stoc, &config.RUN_SEED[0], &cell->omega->vals, &cell->omega->probs, &cell->omega->cnt, config.TOLERANCE);
			for ( m = 0; m < cell->omega->cnt; m++ )
				for ( n = 1; n <= stoc->numOmega; n++ )
					cell->omega->vals[m][n] -= stoc->mean[n-1];
		}

		tic = clock();
		/* Use two-stage algorithm to solve the problem */
		if ( solveBendersCell(stoc, prob, cell) ) {
			errMsg("algorithm", "benders", "failed to solve the cells using MASP algorithm", 0);
			goto TERMINATE;
		}
		cell->time->repTime = ((double) (clock() - tic))/CLOCKS_PER_SEC;

		/* Write solution statistics for optimization process */
		writeStatistic(sFile, iFile, prob, cell);
		writeStatistic(stdout, NULL, prob, cell);

		/* evaluating the optimal solution*/
		if (config.EVAL_FLAG == 1) {
			if ( config.MASTER_TYPE == PROB_QP )
				evaluate(sFile, stoc, prob, cell, cell->incumbX);
			else
				evaluate(sFile, stoc, prob, cell, cell->candidX);
		}
	}

	fclose(sFile); fclose(iFile);
	printf("\nSuccessfully completed the L-shaped method.\n");

	/* free up memory before leaving */
	freeCellType(cell);
	freeProbType(prob, 2);
	mem_free(meanSol);
	return 0;

	TERMINATE:
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	mem_free(meanSol);
	return 1;
}//END algo()

int solveBendersCell(stocType *stoc, probType **prob, cellType *cell) {
	int 	candidCut;
	clock_t	tic, mainTic;

#if defined(SAVE_DUALS)
	if ( duals == NULL ) {
		duals = (dualsType *) mem_malloc(sizeof(dualsType));
		duals->iter = (intvec) arr_alloc(config.MAX_ITER*cell->omega->cnt, int);
		duals->obs = (intvec) arr_alloc(config.MAX_ITER*cell->omega->cnt, int);
		duals->vals = (vector *) arr_alloc(config.MAX_ITER*cell->omega->cnt, vector);
		duals->cnt = 0;
	}
#endif

	mainTic = clock();
	/* Main loop of the algorithm */
	while (cell->k < config.MAX_ITER) {
		tic = clock();

		cell->k++;

#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\nIteration-%d :: Incumbent estimate = %lf; Candidate estimate = %lf.\n", cell->k, cell->incumbEst, cell->candidEst);
#else
		if ( (cell->k-1) % 100 == 0) {
			printf("\nIteration-%4d: ", cell->k);
		}
#endif
		/******* 1a. Optimality tests *******/
		if ( config.MASTER_TYPE == PROB_QP )
			if (optimal(cell))
				break;

		/******* 2. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formOptCut(prob[1], cell, cell->candidX, FALSE)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}

		/******* 4. Check improvement in predicted values at candidate solution *******/
		if ( config.MASTER_TYPE == PROB_QP ) {
			if ( cell->k > 1 ) {
				/* If the incumbent has not changed in the current iteration */
				checkImprovement(prob[0], cell, candidCut);
			}
			else
				cell->incumbEst = vXvSparse(cell->incumbX, prob[0]->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->incumbX, prob[0]->num->cols);
		}
		else {
			cell->incumbEst = vXvSparse(cell->candidX, prob[0]->dBar) + cutHeight(cell->cuts->vals[candidCut], cell->candidX, prob[0]->num->cols);

			if (optimal(cell))
				break;
		}

		/******* 3. Solve the master problem to obtain the new candidate solution */
		if ( solveMaster(prob[0]->num, prob[0]->dBar, cell) ) {
			errMsg("algorithm", "solveCell", "failed to solve master problem", 0);
			return 1;
		}
		cell->time->masterAccumTime += cell->time->masterIter; cell->time->subprobAccumTime += cell->time->subprobIter;
		cell->time->masterIter = cell->time->subprobIter = cell->time->optTestIter = 0.0;
		cell->time->iterTime = ((double) clock() - tic)/CLOCKS_PER_SEC; cell->time->iterAccumTime += cell->time->iterTime;

		if ( cell->k % 10 == 0 )
			printf("Time = %lf", ((double) (clock() - mainTic))/CLOCKS_PER_SEC);
	}

	return 0;
}//END solveCell()

BOOL optimal(cellType *cell) {

	if ( cell->k > config.MIN_ITER ) {
		return cell->optFlag = ((cell->incumbEst - cell->candidEst) < config.EPSILON);
	}

	return FALSE;
}//END optimal()

void writeStatistic(FILE *soln, FILE *incumb, probType **prob, cellType *cell) {

	fprintf(soln, "\n------------------------------------------------------------ Optimization ---------------------------------------------------------\n");
	if ( config.MASTER_TYPE == PROB_QP )
		fprintf(soln, "Algorithm                          : Regularized Benders Decomposition\n");
	else
		fprintf(soln, "Algorithm                          : Benders Decomposition\n");
	fprintf(soln, "Number of iterations               : %d\n", cell->k);
	fprintf(soln, "Lower bound estimate               : %f\n", cell->incumbEst);
	if ( cell->k == config.MAX_ITER)
		fprintf(soln, "Maximum itertions reached with gap.");
	fprintf(soln, "Total time                         : %f\n", cell->time->repTime);
	fprintf(soln, "Total time to solve master         : %f\n", cell->time->masterAccumTime);
	fprintf(soln, "Total time to solve subproblems    : %f\n", cell->time->subprobAccumTime);
	fprintf(soln, "Total time in verifying optimality : %f\n", cell->time->optTestAccumTime);

	if ( incumb != NULL ) {
		printVector(cell->incumbX, prob[0]->num->cols, incumb);
	}

	if ( duals != NULL ) {
		printf("\nNumber of duals discovered = %d", duals->cnt);
	}

}//END WriteStat
