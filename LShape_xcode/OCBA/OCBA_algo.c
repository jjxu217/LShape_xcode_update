/*
 * OCAB_algo.c
 *
 *  Created on: March 16, 2020
 *      Author: Jiajun Xu
 * Institution: University of Southern California
 *
 * Please send your comments or bug report to jiajunx (at) usc (dot) edu
 *
 */
#include <math.h>
#include <stdio.h>
#include "OCBA_algo.h"

extern string outputDir;
extern configType config;

#if defined(SAVE_DUALS)
dualsType *duals = NULL;
#endif

int algo (oneProblem *orig, timeType *tim, stocType *stoc, string probName) {
	probType **prob = NULL;
	cellType *cell = NULL;
	vector 	 meanSol;
    double   std=0, ocba_time=0, naive_time=0, stdev=0, pr, temp;
    batchSummary *batch = NULL;
	int 	 rep, m, n, out_idx=0;
	FILE 	*sFile, *iFile = NULL;
    char results_name[BLOCKSIZE];
    char incumb_name[BLOCKSIZE];
	clock_t	tic;
    
    BOOL distinct;
    double   inverse_appearance[config.NUM_REPS];
    int Delta = 1000;
    int     i, j=0;
    ocbaSummary *ocba = NULL;
    
    
	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell, &batch, &meanSol) )
		goto TERMINATE;

    if ( config.NUM_REPS > 1 )
        ocba  = newOcbaSummary(prob[0]->sp->mac, config.NUM_REPS);
        
	printf("Starting Benders decomposition.\n");
    sprintf(results_name, "results%d.txt", out_idx);
	sFile = openFile(outputDir, results_name, "w");
//    if ( config.MASTER_TYPE == PROB_QP )
    sprintf(incumb_name, "incumb%d.txt", out_idx);
    iFile = openFile(outputDir, incumb_name, "w");
	printDecomposeSummary(sFile, probName, tim, prob);
	printDecomposeSummary(stdout, probName, tim, prob);
    fprintf(sFile, "\n Number of observations in each replications: %d", config.MAX_OBS);

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
//		if (config.EVAL_FLAG == 1) {
//			if ( config.MASTER_TYPE == PROB_QP )
//				evaluate(sFile, stoc, prob, cell, cell->incumbX);
//			else
//				evaluate(sFile, stoc, prob, cell, cell->candidX);
//		}
        /* check is the new solution is distict or not. */
        if ( config.MULTIPLE_REP ) {
            buildCompromise(prob[0], cell, batch);
            
            distinct = TRUE;
            for(i = 0; i < ocba->cnt; i++ ){
                if (equalVector(ocba->incumbX[i], cell->candidX, prob[0]->num->cols, config.TOLERANCE)){
                    distinct = FALSE;
                    break;
                }
            }
            if(distinct == TRUE){
                ocba->objLB[ocba->cnt] = cell->candidEst;
                ocba->incumbX[ocba->cnt] = duplicVector(cell->candidX, prob[0]->num->cols);
                ocba->appearance[ocba->cnt] = 1;
                ocba->cnt++;
            }
            else{
                ocba->appearance[i]++;
                temp = ocba->objLB[i];
                ocba->objLB[i] = (ocba->appearance[i] - 1.0) / ocba->appearance[i] * temp + cell->candidEst / ocba->appearance[i];
            }
        }
	}
    
    if ( config.MULTIPLE_REP ) {
        
        if ( solveCompromise(prob[0], batch)) {
            errMsg("algorithm", "algo", "failed to solve the compromise problem", 0);
            goto TERMINATE;
        }
        
        for (i = 0; i < ocba->cnt; i++)
            inverse_appearance[i] = 1.0 / ocba->appearance[i];
        
        /* Solve the ocba problem. */
        tic = clock();
        
        ocba->idx = solveOCBA(ocba->objLB, inverse_appearance, ocba->cnt, ocba->n, Delta, ocba->an);
        eval_all(sFile, stoc, prob, cell, ocba);
        
        stdev = sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]);
        while(3.92 * stdev > config.EVAL_ERROR * DBL_ABS(ocba->mean[ocba->idx]) || ocba->n[ocba->idx] < config.EVAL_MIN_ITER  ){
            ocba->idx = solveOCBA(ocba->mean, ocba->var, ocba->cnt, ocba->n, Delta, ocba->an);
            eval_all(sFile, stoc, prob, cell, ocba);
        }
        
        ocba_time = ((double) (clock() - tic)) / CLOCKS_PER_SEC;
        batch->time->repTime += ocba_time;
        
        tic = clock();
        for (i = 0; i < ocba->cnt; i++){
            for(j = 0; j < ocba->appearance[i]; j++)
                evaluate(sFile, stoc, prob, cell, ocba->incumbX[i]);
        }
        naive_time = ((double) (clock() - tic)) / CLOCKS_PER_SEC;
        
        fprintf(sFile, "\n====================================================================================================================================\n");
        fprintf(sFile, "\n----------------------------------------- Final solution --------------------------------------\n\n");
        fprintf(sFile, "\n====================================================================================================================================\n");
        fprintf(sFile, "\n----------------------------------------- Final solution --------------------------------------\n\n");
        /* Evaluate the compromise solution */
        fprintf(sFile, "Incumbent solution, non-zero position: ");
        printVectorInSparse(ocba->incumbX[ocba->idx], prob[0]->num->cols, sFile);
        //evaluate(sFile, stoc, prob, cell, ocba->incumbX[ocba->idx]);
        
        
        fprintf(sFile, "Time to solve ocba problem   : %f\n", ocba_time);
        fprintf(sFile, "Total time                         : %f\n", batch->time->repTime);
        fprintf(sFile, "Total time to solve master         : %f\n", batch->time->masterAccumTime);
        fprintf(sFile, "Total time to solve subproblems    : %f\n", batch->time->subprobAccumTime);
        
        fprintf(sFile, "Lower bound estimate               : %f\n", batch->Est);
        std = LowerBoundVariance(batch);
        fprintf(sFile, "Lower bound estimation std               : %f\n", std);
        
//        fprintf(sFile, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
//        fprintf(stdout, "\n------------------------------------------- Average solution ---------------------------------------\n\n");
//        /* Evaluate the average solution */
//        printVectorInSparse(batch->avgX, prob[0]->num->cols, sFile);
//        evaluate(sFile, stoc, prob, cell, batch->avgX);
        
        fprintf(iFile, "\n----------------------------------------- Compromise solution(1-indexed) --------------------------------------\n\n");
        printVectorInSparse(ocba->incumbX[ocba->idx], prob[0]->num->cols, iFile);
        
        /*calculate the   Approximate  Probability  of  Correct  Selection(APCS),
         see: https://math.stackexchange.com/questions/178334/the-probability-of-one-gaussian-larger-than-another */
        pr = 1;
        for (i = 0; i < ocba->cnt; i++){
            if (i != ocba->idx)
                pr -= 0.5 * erfc((ocba->mean[i] - ocba->mean[ocba->idx]) / sqrt(2 * (ocba->var[i] + ocba->var[ocba->idx])));
        }
            
        printf("The Approximate  Probability  of  Correct  Selection(APCS) is %lf", pr);
        
        
        fprintf(sFile, "\n print for latex\n");
        fprintf(sFile, "&%.2f &%.2f &%.2f%% &%.2f &%.2f%% &[%.2f, %.2f], &[%.2f, %.2f] &%.2f", ocba_time, naive_time, 100.0*(naive_time - ocba_time)/naive_time, ocba->objLB[ocba->idx], 100*(1- 1/sqrt(ocba->appearance[ocba->idx])), batch->Est - 1.96*std, batch->Est + 1.96*std, ocba->mean[ocba->idx] - 1.96 * sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]), ocba->mean[ocba->idx] + 1.96 * sqrt(ocba->var[ocba->idx] / ocba->n[ocba->idx]), pr);
        

    }

	fclose(sFile); fclose(iFile);
     
        
//        /*outer loop condition*/
//
//        if (InConvexHull(batch, prob[0]->num->cols) && (std  < config.std_tol)){
//            printf("\nSuccessfully completed the L-shaped method.\n");
//            break;
//        }
//
//        config.MAX_OBS = 2 * config.MAX_OBS;
//        out_idx++;
//    }
	


	/* free up memory before leaving */
	freeCellType(cell);
	freeProbType(prob, 2);
	mem_free(meanSol);
    freeOCBA(ocba, config.NUM_REPS);
	return 0;

	TERMINATE:
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	mem_free(meanSol);
    freeOCBA(ocba, config.NUM_REPS);
	return 1;
}//END algo()

int solveBendersCell(stocType *stoc, probType **prob, cellType *cell) {
    int     candidCut;
    clock_t    tic, mainTic;

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
//        if ( config.MASTER_TYPE == PROB_QP )
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

	if ( cell->RepeatedTime > 0 || cell->k > config.MIN_ITER ) {
        if (config.MASTER_TYPE == PROB_QP || config.reg == 1){
            return cell->optFlag = ((cell->incumbEst - cell->candidEst) < config.EPSILON);
        }
        return TRUE;
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

    //print candidate for a moment here
	if ( incumb != NULL ) {
        printVectorInSparse(cell->candidX, prob[0]->num->cols, incumb);
		//printVector(cell->candidX, prob[0]->num->cols, incumb);
	}

}//END WriteStat

ocbaSummary *newOcbaSummary(int first_stage_cols, int numBatch) {
    ocbaSummary *ocba;
    int i;
    
    ocba = (ocbaSummary *) mem_malloc(sizeof(ocbaSummary));
    
    ocba->cnt = 0;
    ocba->idx = 0;
    
    ocba->ck = (intvec) arr_alloc(numBatch, int);
    ocba->objLB = (vector) arr_alloc(numBatch, double);

    ocba->incumbX = (vector *) arr_alloc(numBatch, vector);
    ocba->appearance = (intvec) arr_alloc(numBatch, int);
    
    ocba->mean = (vector) arr_alloc(numBatch, double);
    ocba->var = (vector) arr_alloc(numBatch, double);
    
    for (i = 0; i < numBatch; i++){
        ocba->incumbX[i] = arr_alloc(first_stage_cols + 1, double);
    }
    
    ocba->n = arr_alloc(numBatch, int);
    ocba->an = arr_alloc(numBatch, int);
    return ocba;
}

void freeOCBA(ocbaSummary *ocba, int numBatch){
    int i;
    
    if ( ocba ) {
        if (ocba->ck) mem_free(ocba->ck);
        if (ocba->objLB) mem_free(ocba->objLB);
        if (ocba->mean) mem_free(ocba->mean);
        if (ocba->var) mem_free(ocba->var);
        if (ocba->appearance) mem_free(ocba->appearance);
        for (i = 0; i < numBatch; i++){
            if (ocba->incumbX[i]) mem_free(ocba->incumbX[i]);
        }
        if (ocba->incumbX) mem_free(ocba->incumbX);
        if(ocba->n) mem_free(ocba->n);
        if(ocba->an) mem_free(ocba->an);
        mem_free(ocba);
    }
}


int solveOCBA(vector s_mean, vector s_var, int nd, intvec n, int add_budget, intvec an){
    /* This subroutine determines how many additional runs each design will should have for next iteration of simulation.
    s_mean[i]: sample mean of design i, i=0,1,..,nd-1
    s_var[i]: sample variance of design i, i=0,1,..,nd-1
    nd: the number of designs
    n[i]: number of simulation replications of design i, i=0,1,..,nd-1
    add_budget: the additional simulation budget
    an[i]: additional number of simulation replications assigned to design i, i=0,1,..,nd-1 */
    
    int i;
    int b, s;
    int t_budget, t1_budget;
    int morerun[nd], more_alloc; /* 1:Yes; 0:No */
    double t_s_mean[nd];
    double ratio[nd]; /* Ni/Ns */
    double ratio_s, temp;
    
    for(i = 0; i < nd; i++)
        t_s_mean[i] = s_mean[i];
    
    /*t_budget record total budget*/
    t_budget = add_budget;
    for(i = 0; i < nd; i++)
        t_budget += n[i];
    b = best(t_s_mean, nd);
    s = second_best(t_s_mean, nd, b);

    /* calculate ratio of Ni/Ns*/
    ratio[s] = 1.0;
    for(i = 0; i < nd; i++)
        if(i != s && i != b){
            temp = (t_s_mean[b] - t_s_mean[s]) / (t_s_mean[b] - t_s_mean[i]);
            ratio[i] = temp * temp * s_var[i] / s_var[s];
        }
    
    /* calculate Nb */
    temp = 0;
    for(i = 0; i < nd; i++)
        if(i != b)
            temp += (ratio[i] * ratio[i] / s_var[i]);
    ratio[b] = sqrt(s_var[b] * temp);
    
    
    for(i = 0; i < nd; i++)
        morerun[i] = 1;
    t1_budget = t_budget;
    
    do{
        more_alloc = 0;
        ratio_s = 0.0;
        for(i = 0; i < nd; i++)
            if(morerun[i])
                ratio_s += ratio[i];
        for(i = 0; i < nd; i++)
            if(morerun[i]) {
                an[i] = (int)(t1_budget / ratio_s * ratio[i]);
                
                /* disable those design which have been run too much */
                if(an[i] <= n[i]){
                    an[i] = n[i];
                    morerun[i] = 0;
                    more_alloc = 1;
                }
            }
        if(more_alloc) {
            t1_budget = t_budget;
            for(i = 0; i < nd; i++)
                if(!morerun[i])
                    t1_budget -= an[i];
        }
    } while(more_alloc); /* end of WHILE */
    
    /* calculate the difference */
    t1_budget = 0;
    for(i = 0; i < nd; i++)
        t1_budget += an[i];
    an[b] += (t_budget - t1_budget); /* give the difference to design b */
    for(i = 0; i < nd; i++)
        an[i] -= n[i];
    return b;
}



int eval_all(FILE *soln, stocType *stoc, probType **prob, cellType *cell, ocbaSummary *ocba) {
    vector     observ, rhs, costTemp, cost;
    intvec    objxIdx;
    double     obj, mean, variance, stdev, temp, pre_mean;
    int        cnt, m, status, i;

    if ( !(observ = (vector) arr_alloc(stoc->numOmega + 1, double)) )
        errMsg("allocation", "evaluateOpt", "observ", 0);
    
    /* right-hand side */
    if (!(rhs =(vector) arr_alloc(prob[1]->num->rows+1, double)))
        errMsg("Allocation", "evaluate", "rhs",0);

    /* cost coefficients */
    if ( !(cost = (vector) arr_alloc(prob[1]->num->cols+1, double)) )
        errMsg("allocation", "evaluate", "cost", 0);
    if ( !(objxIdx = (intvec) arr_alloc(prob[1]->num->cols+1, int)) )
        errMsg("allocation", "evaluate", "objxIdx", 0);
    costTemp = expandVector(prob[1]->dBar->val, prob[1]->dBar->col, prob[1]->dBar->cnt, prob[1]->num->cols);
    for (m = 1; m <= prob[1]->num->rvdOmCnt; m++ ) {
        objxIdx[m] = prob[1]->coord->rvdOmCols[m] - 1;
        cost[m] = costTemp[objxIdx[m]+1];
    }
    mem_free(costTemp);
    

    for (i = 0; i < ocba->cnt; i++){
    
        printf("\nStarting evaluate solution %d with obs %d.\n", i, ocba->an[i]);

        /* initialize parameters used for evaluations */
        cnt = 0; mean = 0.0; variance = 0.0; pre_mean = 0; //stdev = INFBOUND;


        /* change the right hand side with the solution */
        chgRHSwSoln(prob[1]->bBar, prob[1]->Cbar, rhs, ocba->incumbX[i]);

        while (cnt < ocba->an[i] ) {
            /* use the stoc file to generate observations */
            generateOmega(stoc, observ, config.TOLERANCE, &config.EVAL_SEED[0]);

            for ( m = 0; m < stoc->numOmega; m++ )
                observ[m] -= stoc->mean[m];

            /* Change right-hand side with random observation */
            if ( chgRHSwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, observ-1, rhs, ocba->incumbX[i]) ) {
                errMsg("algorithm", "evaluate", "failed to change right-hand side with random observations",0);
                return 1;
            }

            /* Change cost coefficients with random observations */
            if ( prob[1]->num->rvdOmCnt > 0 ) {
                if ( chgObjxwObserv(cell->subprob->lp, prob[1]->num, prob[1]->coord, cost, objxIdx, observ-1) ) {
                    errMsg("algorithm", "evaluate","failed to change cost coefficients with random observations", 0);
                    return 1;
                }
            }

            if ( solveProblem(cell->subprob->lp, cell->subprob->name, cell->subprob->type, &status) ) {
                if ( status == STAT_INFEASIBLE ) {
                    /* subproblem is infeasible */
                    printf("Warning:: Subproblem is infeasible: need to create feasibility cut.\n");
                    return 1;
                }
                else {
                    errMsg("algorithm", "evaluateOpt", "failed to solve subproblem in solver", 0);
                    return 1;
                }
            }

            /* use subproblem objective and compute evaluation statistics */
            obj = getObjective(cell->subprob->lp, PROB_LP);

    #if defined(ALGO_CHECK)
            writeProblem(cell->subprob->lp, "evalSubprob.lp");
            printf("Evaluation objective function = %lf.\n", obj);
    #endif


            if ( cnt == 0 )
                mean = obj;
            else {
                temp = mean;
                mean = mean + (obj - mean) / (double) (cnt + 1);
                variance  = cnt / (cnt + 1) * variance
                + cnt * (mean - temp) * (mean - temp);
                //stdev = sqrt(variance/ (double) cnt);
            }
            cnt++;

            /* Print the results every once in a while for long runs */
            if (!(cnt % 100)) {
                printf(".");
                fflush(stdout);
            }
            if (!(cnt % 10000))
                printf("\nObs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", cnt, mean, 3.92 * stdev / mean,  mean - 1.96 * stdev, mean + 1.96 * stdev);
            
    
        }//END while loop
        mean += vXvSparse(ocba->incumbX[i], prob[0]->dBar);
        
        /* New mean and variance
         see: https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-or-more-groups-given-known-group-varianc
        \sigma^2_{1:m+n} = \frac{n(\sigma^2_{1:n} + \mu_{1:n}^2) + m(\sigma^2_{1+n:m+n} + \mu_{1+n:m+n}^2)}{m+n} - \mu^2_{1:m+n}.*/
        
        pre_mean = ocba->mean[i];
        ocba->mean[i] = (ocba->n[i] * pre_mean + ocba->an[i] * mean) / (ocba->n[i] + ocba->an[i]);
        
        ocba->var[i] = (double) ocba->n[i] / (ocba->n[i] + ocba->an[i]) * (ocba->var[i] + pre_mean * pre_mean) + (double) ocba->an[i] / (ocba->n[i] + ocba->an[i]) * (variance + mean * mean) - ocba->mean[i] * ocba->mean[i];
        
        ocba->n[i] += ocba->an[i];
        ocba->an[i] = 0;
        
        stdev = sqrt(ocba->var[i] / ocba->n[i]);
        printf("\n mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n", ocba->mean[i], 3.92 * stdev / ocba->mean[i],  ocba->mean[i] - 1.96 * stdev, ocba->mean[i] + 1.96 * stdev);
        
        

//        printf("\n\nEvaluation complete. Final evaluation results :: \n");
//        printf("Upper bound estimate               : %lf\n", mean);
//        printf("Error in estimation                : %lf\n", 3.92 * stdev / mean);
//        printf("Confidence interval at 95%%         : [%lf, %lf]\n", mean - 1.96 * stdev, mean + 1.96 * stdev);
//        printf("Number of observations             : %d\n", cnt);
//
//        if ( soln != NULL ) {
//            /* Write the evaluation results to the summary file */
//            fprintf(soln, "------------------------------------------------------------- Evaluation ----------------------------------------------------------\n");
//            fprintf(soln, "Upper bound estimate           : %lf\n", mean);
//            fprintf(soln, "Error in estimation            : %lf\n", 3.92 * stdev / mean);
//            fprintf(soln, "Confidence interval at 95%%     : [%lf, %lf]\n", mean - 1.96 * stdev, mean + 1.96 * stdev);
//            fprintf(soln, "Number of observations         : %d\n", cnt);
//        }

    }
    mem_free(observ); mem_free(rhs); mem_free(objxIdx); mem_free(cost);
    return 0;

}//END evaluate()

int best(vector t_s_mean, int nd){
    /*This function determines the best design based on current simulation results
     t_s_mean[i]: temporary array for sample mean of design i, i=0,1,..,ND-1
     nd: the number of designs */
    int i, min_index=0;
    for (i = 0; i < nd; i++){
        if(t_s_mean[i] < t_s_mean[min_index])
            min_index = i;
    }
    return min_index;
}

int second_best(vector t_s_mean, int nd, int b){
    int i, second_index;
    if (b==0)
        second_index = 1;
    else
        second_index = 0;
    for(i = 0;i < nd; i++){
        if(t_s_mean[i] < t_s_mean[second_index] && i != b)
        {
            second_index = i;
        }
    }
    return second_index;
}
