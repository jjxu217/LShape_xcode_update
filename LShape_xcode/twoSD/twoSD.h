/*
 * twoSD.h
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef TWOSD_H_
#define TWOSD_H_

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"
#include "stoc.h"

#define TRIVIAL 0
#define NONTRIVIAL 1

#undef ALGO_CHECK
#undef BATCH_CHECK

/* A data structure which holds on the configuration information about the algorithm. Most of these configuration parameters are read from a
-configuration file. These elements, once set during initialization, are not modified during the course of the algorithm. */
typedef struct{
	int		NUM_REPS;			/* Maximum number of replications that can be carried out. */
	long long *RUN_SEED;		/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTER_TYPE;		/* type of master problem */
	int		TAU;				/* Frequency at which the incumbent is updated */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double	EPSILON;			/* Optimality gap */

	int		EVAL_FLAG;
	int		NUM_EVALS;
	long long *EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;

	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double	R1;
	double	R2;
	double	R3;
	int		DUAL_STABILITY;		/* Determine whether dual stability is to be checked. */
	int		PI_EVAL_START;		/* The minimum number of samples before dual stability test is conducted. */
	int		PI_CYCLE;			/* Frequency of updating the dual stability ratio */
	int		BOOTSTRAP_REP;		/* Number of boot-strap replications in full optimality test */
	double	PERCENT_PASS;		/* percentage of bootstrap replications need to be satisfied */
	int		SCAN_LEN;			/* window size over which the stability of dual vertex set is measured.*/
	double  PRE_EPSILON;		/* gap used for preliminary optimality test */

	int 	MULTIPLE_REP;		/* When multiple replications are needed, set this to (M), else (0) */
}configType;

typedef struct {
	double  alpha;                  /* scalar value for the righ-hand side */
	vector  beta;                   /* coefficients of the master problems's primal variables */
	int 	numSamples;				/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut is based (this is also the length of istar) */
	intvec	iStar;					/* indices of maximal pi for each distint observation */
	BOOL	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	slackCnt;				/* number of times a cut has been slack, used in deciding when the cut needs to be dropped */
	int 	rowNum;					/* row number for master problem in solver */
	string	name;
}oneCut;

typedef struct {
	int     cnt;                    /* number of cuts */
	oneCut  **vals;
}cutsType;

typedef struct {
	double	repTime;
	double 	iterTime;
	double 	masterIter;
	double 	subprobIter;
	double 	optTestIter;
	double 	argmaxIter;
	double 	iterAccumTime;
	double 	masterAccumTime;
	double 	subprobAccumTime;
	double 	optTestAccumTime;
	double 	argmaxAccumTime;
}runTime;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */
    double		lb;					/* lower bound on cell objective function */
    int			lbType;				/* type of lower bound being used TRIVIAL if 0, else NONTRIVIAL */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	vector      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	vector      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */

	vector 		piM;				/* master dual information */
	vector      djM;                /* master reduced cost vector */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fcuts;             /* feasibility cuts */
    cutsType 	*fcutsPool;			/* Pool of feasibility cuts */
    int			fUpdt[2];			/* coordinate in delta structure for which the updates have been carried out */

	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */
	basisType	*basis;				/* hold unique basis identified */

    BOOL        optFlag;
	vector      pi_ratio;
    BOOL        dualStableFlag; 	/* indicates if dual variables are stable */

    BOOL 		optMode;			/* When false, the algorithm tries to resolve infeasibility */

    BOOL		spFeasFlag;			/* Indicates whether the subproblem is feasible */
	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasible */

	runTime		time;				/* Run time structure */
}cellType;

typedef struct {
	oneProblem	*sp;				/* compromise problem */
	int 		cnt;				/* number of replications */
	intvec 		ck;					/* number of iterations for each replication */
	vector		objLB;				/* replication lower bound */
	vector		objUB;				/* replication upper bound, if batch solution is evaluated */
	double		objComp;			/* optimal value of compromise problem */
	double		quadScalar;			/* average proximal terms */
	vector		*incumbX;			/* batch incumbent solution */
	vector		compromiseX;		/* compromise solution */
	vector		avgX;				/* average solution across batches */
}batchSummary;

/* twoSD.c */
void parseCmdLine(int argc, char *argv[], string probName, string inputDir);
int readConfig();

/* algo.c */
int algo(oneProblem *orig, timeType *tim, stocType *stoc, string inputDir, string probName);
int solveCell(stocType *stoc, probType **prob, cellType *cell);
void writeOptimizationSummary(FILE *soln, FILE *incumb, probType **prob, cellType *cell, BOOL header);
void cleanupAlgo(probType **prob, cellType *cell, int T);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell, batchSummary **batch, vector *meanSol);
cellType *newCell(stocType *stoc, probType **prob, vector xk);
int cleanCellType(cellType *cell, probType *prob, vector xk);
void freeCellType(cellType *cell);

/* master.c */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, double lb);
int addCut2Master(oneProblem *master, oneCut *cut, vector vectX, int lenX);
int constructQP(probType *prob, cellType *cell, vector incumbX, double quadScalar);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell, vector xk);
int changeQPbds(LPptr lp, int numCols, vector bdl, vector bdu, vector xk, int offset);
oneProblem *newMaster(oneProblem *orig, double lb);

/* cuts.c */
int formSDCut(probType **prob, cellType *cell, vector Xvect, int omegaIdx, BOOL *newOmegaFlag, double lb);
oneCut *SDCut(numType *num, coordType *coord, basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, vector Xvect, int numSamples,
		BOOL *dualStableFlag, vector pi_ratio, double lb);
oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
int reduceCuts(cellType *cell, vector candidX, vector pi, int betaLen, double lb);
int dropCut(cellType *cell, int cutIdx);
double calcVariance(double *x, double *mean_value, double *stdev_value, int batch_size);
void printCut(cutsType *cuts, numType *num, int idx);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts, BOOL partial);
double calc_var(double *x, double *mean_value, double *stdev_value, int batch_size);

/* soln.c */
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell, double candidEst);
double maxCutHeight(cutsType *cuts, int currIter, vector xk, int betaLen, double lb);
double cutHeight(oneCut *cut, int currIter, vector xk, int betaLen, double lb);

/* optimal.c */
BOOL optimal(probType **prob, cellType *cell);
BOOL preTest(cellType *cell);
BOOL fullTest(probType **prob, cellType *cell);
cutsType *chooseCuts(cutsType *cuts, vector pi, int lenX);
void reformCuts(basisType *basis, sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord,
		cutsType *gCuts, int *observ, int k, int lbType, int lb, int lenX);
double calcBootstrpLB(probType *prob, vector incumbX, vector piM, vector djM, int currIter, double quadScalar, cutsType *cuts);
void empiricalDistribution(omegaType *omega, int *cdf);
void resampleOmega(intvec cdf, intvec observ, int numSamples);

/* compromise.c */
int buildCompromise(probType *prob, cellType *cell, batchSummary *batch);
int solveCompromise(probType *prob, batchSummary *batch);
int addBatchEquality (probType *prob, batchSummary *batch);
batchSummary *newBatchSummary(probType *prob, int numBatch);
void freeBatchType(batchSummary *batch);

/* evaluate.c */
int evaluate(FILE *soln, stocType *stoc, probType **prob, oneProblem *subprob, vector Xvect);
void writeEvaluationSummary(FILE *soln, double mean, double stdev, int cnt);

#endif /* TWOSD_H_ */
