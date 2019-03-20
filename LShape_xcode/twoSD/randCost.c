/*
 * randCost.c
 *
 *  Created on: Jan 10, 2018
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"
#include "stoc.h"

extern string outputDir;

void calcBasis(LPptr lp, numType *num, coordType *coord, sparseVector *dBar, oneBasis *B, int basisDim) {
	vector	basicCost, costVector, tempPsiRow, basisVal;
	intvec 	basisHead, phiHead;
	int		i, j;

	/* Allocate memory for the basis header. */
	basisHead = (intvec) arr_alloc(basisDim+1, int);
	phiHead = (intvec) arr_alloc(num->rvdOmCnt+1, int);
	B->gBar = (vector) arr_alloc(num->cols+1, double);

	/* Compute the phi matrix associated with the current basis. We begin by first identifying the basis header. A negative
	 * value in basis header indicates a slack row. */
	getBasisHead(lp, basisHead+1, NULL);

	/* Compute the phi matrix header and extract the basis (of the dual) inverse matrix rows corresponding to the header. */
	for ( i = 1; i <= num->rvdOmCnt; i++ ) {		/* Loop through all the columns with random cost coefficients to see if any of them are basic */
		if ( (phiHead[i] = isElementIntvec(basisHead, basisDim, (coord->rvdOmCols[i]-1) )) > 0 ) {
			/* _idx_ > 0 gives the index of column with random cost coefficient in the basis head */
			if ( B->phiLength == 0 ) {
				if ( !(B->phi = (vector *) arr_alloc(num->rvdOmCnt, vector)) )
					errMsg("allocation", "newBasis", "B->phi", 0);
				if ( !(B->omegaIdx = (intvec) arr_alloc(num->rvdOmCnt+1, int)) )
					errMsg("allocation", "newBasis", "B->lambdaIdx", 0);
			}
			if ( !(B->phi[B->phiLength] = (vector) arr_alloc(basisDim+1, double)) )
				errMsg("allocation", "calcBasis", "B->phi[i]", 0);
			getBasisInvRow(lp, phiHead[i]-1, B->phi[B->phiLength]+1);
			B->omegaIdx[++B->phiLength] = i;
		}
	}

	if ( B->phiLength > 0 ) {
		/* Reallocate memory to elements already assigned and initialize for the remainder of the elements. */
		B->phi = (vector *) mem_realloc(B->phi, B->phiLength*sizeof(vector));

		B->omegaIdx = (intvec) mem_realloc(B->omegaIdx, (B->phiLength+1)*sizeof(int));
		B->sigmaIdx = (intvec) mem_realloc(B->sigmaIdx, (B->phiLength+1)*sizeof(int));

		if ( !(B->psi = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
			errMsg("allocation", "newBasis", "B->psi", 0);
		B->psi->val = (vector) arr_alloc(num->cols*B->phiLength+1, double);
		B->psi->col = (intvec) arr_alloc(num->cols*B->phiLength+1, int);
		B->psi->row = (intvec) arr_alloc(num->cols*B->phiLength+1, int);
		B->psi->cnt = 0;
	}

	/* Extract the basic variable cost vector and psi matrix (the tableau entries) */
	if ( !(basicCost = (vector) arr_alloc(basisDim+1, double)) )
		errMsg("allocation", "newBasis", "basicCost", 0);
	costVector = expandVector(dBar->val, dBar->col, dBar->cnt, num->cols);
	for ( i = 1; i <= basisDim; i++ ) {
		if ( basisHead[i] >= 0 )
			basicCost[i] = costVector[basisHead[i]+1];
	}

	if ( !(tempPsiRow = (vector) arr_alloc(num->rows+1, double)) )
		errMsg("allocation", "newBasis", "tempPsiRow", 0);

	for ( i = 1; i <= num->cols; i++ ) {
		getBasisInvACol(lp, i-1, tempPsiRow+1);

		B->gBar[i] = costVector[i] - vXv(tempPsiRow, basicCost, NULL, num->rows);

		for ( j = 1; j <= B->phiLength; j++ ) {
			B->psi->row[B->psi->cnt+1] = i;
			B->psi->col[B->psi->cnt+1] = B->omegaIdx[j];
			B->psi->val[B->psi->cnt+1] = tempPsiRow[phiHead[B->omegaIdx[j]]];
			B->psi->cnt++;
		}
	}

#if defined (STOCH_CHECK)
	printf("New basis identified     :: \n");
	printf("\tNumber of basic columns with random cost coefficients      = %d\n", B->phiLength);
	if ( B->phiLength > 0 ) {
		printf("\tIndex in observation vector corresponding to basic columns = "); printIntvec(B->omegaIdx, B->phiLength, NULL);
		printf("\tPhi = ");
		for (int i = 0; i < B->phiLength; i++ ) {
			printf("\t\t"); printVector(B->phi[i], num->rows, NULL);
		}
	}
	else {
		printf("\tBasic stochastic columns                                   = NULL\n");
		printf("\tPhi = NULL\n");
	}
	printf("Deterministic component of reduced cost  = ");
	printSparseVector(B->gBar+1, basisHead, num->rows);
#endif

#if defined(BASIS_CHECK)
	FILE *bFile; vector temp;
	bFile = openFile(outputDir, "basis.txt", "w");
	temp = (vector) arr_alloc(num->rows+1, double);
	printIntvec(basisHead, num->rows, bFile);
	printVector(basicCost, num->rows, bFile);
	for ( i = 0; i < num->rows; i++ ) {
		getBasisInvRow(lp, i, temp+1);
		printVector(temp, num->rows, bFile);
	}
	mem_free(temp); fclose(bFile);
#endif

	mem_free(basisHead); mem_free(phiHead); mem_free(costVector); mem_free(basicCost); mem_free(tempPsiRow);

}//END calcBasis()

oneBasis *newBasis(LPptr lp, int numCols, int numRows, int currentIter, BOOL subFeasFlag) {
	oneBasis *B;
	intvec 	cstat, rstat;

	/* allocate memory to elements of the basis structure */
	if ( !(B = (oneBasis *) mem_malloc(sizeof(oneBasis))))
		errMsg("allocation", "newBasis", "B", 0);
	B->sigmaIdx = (intvec) arr_alloc(1, int);
	B->ck    	 = currentIter;
	B->weight 	 = 1;
	B->phiLength = 0;
	B->phi 		 = NULL; B->omegaIdx = NULL; B->gBar = NULL; B->piDet = NULL; B->psi = NULL;
	B->feasFlag = subFeasFlag;

	/* Allocate memory. */
	if ( !(cstat = (intvec) arr_alloc( numCols+1, int)))
		errMsg("allocation", "stochasticUpdates", "cstat", 0);
	if ( !(rstat = (intvec) arr_alloc( numRows+1, int)))
		errMsg("allocation", "stochasticUpdates", "rstat", 0);

	/* Obtain the status of columns and rows in the basis. */
	if ( getBasis(lp, cstat, rstat) ) {
		errMsg("algorithm", "newBasis", "failed to get the basis column and row status", 0);
		return NULL;
	}

#if defined(BASIS_CHECK)
	FILE *basisFile;

	basisFile = openFile(outputDir, "cstat.txt", "a");
	printIntvec(cstat, numCols, basisFile);
	fclose(basisFile);

	basisFile = openFile(outputDir, "rstat.txt", "a");
	printIntvec(cstat, numRows, basisFile);
	fclose(basisFile);

#endif

	if ( computeMU(lp, cstat,  numCols, &B->mubBar) ) {
		errMsg("algorithm", "newBasis", "failed to compute mubBar for subproblem", 0);
		return NULL;
	}

	if ( subFeasFlag ) {
		/* encode the row and column status */
		B->cCode = encodeIntvec(cstat, numCols, WORDLENGTH, 3);
		B->rCode = encodeIntvec(rstat, numRows, WORDLENGTH, 3);
	}
	else {
		B->cCode = B->rCode = NULL;
	}

	mem_free(cstat); mem_free(rstat);
	return B;
}//END newBasis()

int decomposeDualSolution(LPptr lp, oneBasis *B, vector omegaVals, int numRows) {
	int n, i;

	if ( !(B->piDet = (vector) arr_alloc(numRows+1, double)) )
		errMsg("allocation", "decomposeDualSolution", "piS", 0);

	/* Record the dual and reduced cost on bounds. */
	if ( getDual(lp, B->piDet, numRows) ) {
		errMsg("algorithm", "stochasticUpdates", "failed to get the dual", 0);
		return 1;
	}

	for ( n = 0; n < B->phiLength; n++ )
		for ( i = 1; i <= numRows; i++ )
			B->piDet[i] -= B->phi[n][i]*omegaVals[B->omegaIdx[n+1]];

	return 0;
}//END decomposeDualSolution()

/* This subroutine determines the feasibility of a basis. */
BOOL checkBasisFeasibility(oneBasis *B, sparseVector dOmega, string senx, int numCols, int numRows, double TOLERANCE) {
	vector 	reducedCost, theta = NULL;
	intvec	cstat;
	int 	n, c;

	/* If there are no random cost coefficients, then the basis is always feasible. */
	if ( dOmega.cnt > 0 ) {
		/* Reconstruct the dual solution from the basis inverse matrix, and make sure the signs correspond to those expected from
		 * inequality constraints in the primal linear program */
		if ( B->phiLength > 0 ) {
			/* *theta* holds the stochastic component of the dual solution */
			theta = (vector) arr_alloc(numRows+1, double);
			for ( c = 1; c <= numRows; c++ ) {
				for ( n = 0; n < B->phiLength; n++ )
					theta[c] += B->phi[n][c]*dOmega.val[B->omegaIdx[n+1]];
				if ( (B->piDet[c] + theta[c] < -TOLERANCE && senx[c-1] == 'G') || (B->piDet[c] + theta[c] > TOLERANCE && senx[c-1] == 'L')) {
					mem_free(theta);
					return FALSE;
				}
			}
		}

#if defined(BASIS_CHECK)
		printf("Re-constructed dual solution: \n");
		printVector(B->piDet, numRows, NULL);
		if ( B->phiLength > 0 )
			printVector(theta, numRows, NULL);
#endif

		mem_free(theta);

		/* Compute the reduced cost */
		if ( !(reducedCost = (vector) arr_alloc(numCols+1, double)) )
			errMsg("allocation", "calcDelta", "costVector", 0);

		copyVector(B->gBar, reducedCost, numCols, TRUE);
		addVectors(reducedCost, dOmega.val, dOmega.col, dOmega.cnt);
		if ( B->phiLength > 0 ) {
			MSparsexvSub(B->psi, dOmega.val, reducedCost);
		}

		/* Decode the basis column status to determine how the reduced cost must be evaluated */
		cstat = decodeIntvec(B->cCode, numCols, WORDLENGTH, 3);
		c = 1;
		while ( c <= numCols ) {
			if ( reducedCost[c] < -TOLERANCE && cstat[c] != AT_UPPER ) {
				mem_free(reducedCost); mem_free(cstat);
				return FALSE;
			}
			c++;
		}

		mem_free(reducedCost); mem_free(cstat);
	}

	return TRUE;
}//END checkBasisFeasibility()

/* This function allocates a new basisType data structure which holds all the unique basis discovered by the algorithm. It returns a pointer to the
 * structure. */
basisType *newBasisType(int numIter, int numCols, int numRows, int wordLength) {
	basisType *basis;
	int numBits = 2; 		/* The column or row status in a basis is indicated by an integer 0,1,2, or 3. We need 2 bits to encode this information. */

	if ( !(basis = (basisType *) mem_malloc(sizeof(basisType))))
		errMsg("allocation", "newBasisType", "basis", 0);
	if ( !(basis->vals = (oneBasis **) arr_alloc(numIter, oneBasis *)))
		errMsg("allocation", "newBasisType", "basis->vals", 0);
	if ( !(basis->obsFeasible = (BOOL **) arr_alloc(numIter, BOOL *)))
		errMsg("allocation", "newBasisType", "basis->obsFeasible", 0);
	basis->cnt = 0;
	basis->basisDim = min(numRows, numCols);
	basis->cCodeLen = ceil((double) numBits*numCols/ (double) wordLength);
	basis->rCodeLen = ceil((double) numBits*numRows/(double) wordLength);

	return basis;
}//END newBasis()

void freeOneBasis(oneBasis *B) {
	int n;

	if ( B ) {
		if (B->cCode) mem_free(B->cCode);
		if (B->rCode) mem_free(B->rCode);
		if (B->sigmaIdx) mem_free(B->sigmaIdx);
		if (B->omegaIdx) mem_free(B->omegaIdx);
		if (B->gBar) mem_free(B->gBar);
		if (B->piDet) mem_free(B->piDet);
		if ( B->phi) {
			for ( n = 0; n < B->phiLength; n++ )
				if (B->phi[n]) mem_free(B->phi[n]);
			mem_free(B->phi);
		}
		if (B->psi) freeSparseMatrix(B->psi);
		mem_free(B);
	}

}//END freeOneBasis

void freeBasisType(basisType *basis, BOOL partial) {
	int n;

	if ( basis ) {
		if ( basis->vals ) {
			for ( n = 0; n < basis->cnt; n++ ) {
				freeOneBasis(basis->vals[n]);
				mem_free(basis->obsFeasible[n]);
			}
			if ( partial ) {
				basis->cnt = 0;
				return;
			}
			mem_free(basis->vals);
			mem_free(basis->obsFeasible);
		}
		mem_free(basis);
	}

}//END freeBasisType
