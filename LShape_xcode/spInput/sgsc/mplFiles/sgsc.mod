# Set Definition 
set K;    # set of facilities (including plants and distribution centers) indexed by k
set K_P;  # set of plants indexed by k
set K_DC; # set of distribution centers indexed by k
set R;    # set of customers indexed by r
set J;    # set of products indexed by j
set M;    # set of transportation mode indexed by m				
set T;    # set of time periods indexed by t
set P2DC_LINKS within {K, K};  # set of transportation links between plants and distribution centers
set DC2R_LINKS within {K, R} ;  # set of transportation links between distribution centers and customers
set P2R_LINKS within {K, R} ; # set of transportation links between plants and customers

set PDC2R_LINKS := DC2R_LINKS union P2R_LINKS ;

# Start working here

# Parameter Definition
# n indicates the total number of time periods 1..n
param n integer > 0;

# demand of product j in customer r at time period t
param d {r in R, j in J, t in T};

# initial inventory level of product j at facility k
param I0 {k in K, j in J};

# freight rate of product j from facility k1 to k2 with mode m at time period t
param gammaK {k1 in K, k2 in K, j in J, m in M, t in T : (k1, k2) in P2DC_LINKS};

# freight rate of product j from facility k to customer r with mode m at time period t
param gammaR {k in K, r in R, j in J, m in M, t in T : (k, r) in PDC2R_LINKS};

# unit inventory cost of product j in facility k at time period t
param h {k in K, j in J, t in T};

# unit throughput cost of product j in facility k at time period t
param delta {k in K, j in J, t in T};

# unit penalty cost of product j for lost unmeet demand in cus- tomer r at time period t
param eta {r in R, j in J, t in T};

# capacity of plant k for product j, k in K_P
param Q {k in K_P, j in J};

# minimum inventory of product j at facility k at time period t
param Im {k in K, j in J, t in T};

# shipping time of product j from facility k1 to facility k2 with mode m
param lambdaK {k1 in K, k2 in K, j in J, m in M : (k1, k2) in P2DC_LINKS};

# shipping time of product j from facility k to customer r with mode m
param lambdaR {k in K, r in R, j in J, m in M : (k, r) in PDC2R_LINKS};


# Variable Definition

# inter-facility freight of product j from facility k1 to k2 with mode m at time period t
var F {k1 in K, k2 in K, j in J, m in M, t in T : (k1, k2) in P2DC_LINKS} >= 0;

# inventory level of product j at facility k at time period t
var I {k in K, j in J, t in T} >= 0;

# facility-customer freight of product j from facility k to customer r with mode m at time period t
var S {k in K, r in R, j in J, m in M, t in T : (k, r) in PDC2R_LINKS} >= 0 ;

# production amount of product j at plant k at time period t
var W {k in K_P, j in J, t in T} >= 0;

# unmeet demand of product j in customer r at time period t
var SF {r in R, j in J, t in T} >= 0;

# Objective Definition
minimize totalCost:
   sum {k in K, j in J, t in T : t = 1} h[k, j, t] * I[k, j, t] + sum {k1 in K, k2 in K, j in J, m in M, t in T : t = 1 and (k1, k2) in P2DC_LINKS} gammaK[k1, k2, j, m, t] * F[k1, k2, j, m, t]
   + sum {k in K, r in R, j in J, m in M, t in T : t = 1 and (k, r) in  PDC2R_LINKS} gammaR[k, r, j, m, t] * S[k, r, j, m, t]
   + sum {k1 in K, k2 in K, j in J, m in M, t in T: t = 1 and (k1, k2) in P2DC_LINKS} delta[k1, j, t] * F[k1, k2, j, m, t]
   + sum {k in K, r in R, j in J, m in M, t in T: t = 1 and (k, r) in DC2R_LINKS} delta[k, j, t] * S[k, r, j, m, t]
   + sum {k in K, j in J, t in 2..n} h[k, j, t] * I[k, j, t] + sum {k1 in K, k2 in K, j in J, m in M, t in 2..n : (k1, k2) in P2DC_LINKS} gammaK[k1, k2, j, m, t] * F[k1, k2, j, m, t]
   + sum {k in K, r in R, j in J, m in M, t in 2..n: (k, r) in PDC2R_LINKS} gammaR[k, r, j, m, t] * S[k, r, j, m, t]
   + sum {k1 in K, k2 in K, j in J, m in M, t in 2..n : (k1, k2) in P2DC_LINKS} delta[k1, j, t] * F[k1, k2, j, m, t]
   + sum{k in K, r in R, j in J, m in M, t in 2..n : (k, r) in DC2R_LINKS} delta[k, j, t] * S[k, r, j, m, t] + sum {r in R, j in J, t in T} eta[r, j, t] * SF[r, j, t] ;

# Constraints Definition

# First stage mass balance for plants
subject to massBalanceP1 {j in J, k1 in K_P, t in T : t = 1} :
   sum {k2 in K, m in M : (k1, k2) in P2DC_LINKS} F [k1, k2, j, m, t] + sum {r in R, m in M : (k1, r) in P2R_LINKS} S [k1, r, j , m, t] = I0 [k1, j] - I [k1, j, t] + W [k1, j, t] ;

# First stage mass balance for distribution centers
subject to massBalanceDC1 {j in J, k1 in K_DC, t in T : t = 1} :
   sum {k2 in K, m in M : (k1, k2) in P2DC_LINKS} F [k1, k2, j, m, t] + sum {r in R, m in M : (k1, r) in DC2R_LINKS} S [k1, r, j, m, t] = I0 [k1, j] - I [k1, j, t] + sum {k2 in K, m in M : (k2, k1) in P2DC_LINKS and t > lambdaK[k2, k1, j, m]} F [k2, k1, j, m, t - lambdaK[k2, k1, j, m]];

# First stage capacity constraints
subject to capacity1 {k in K_P, j in J, t in T : t = 1}:
   W[k, j, t] <= Q[k, j];

# First stage Minimum inventory constraints
subject to minimumInv1 {k in K, j in J, t in T : t = 1}:
   I[k, j, t] >= Im[k, j, t];

# Second stage mass balance for plants
subject to massBalanceP2 {j in J, k1 in K_P, t in 2..n} :
   sum {k2 in K, m in M : (k1, k2) in P2DC_LINKS} F [k1, k2, j, m, t] + sum {r in R, m in M : (k1, r) in PDC2R_LINKS} S [k1, r, j , m, t] = I [k1, j, t-1] - I [k1, j, t] + W [k1, j, t] + sum {k2 in K, m in M : (k2, k1) in P2DC_LINKS} F [k2, k1, j, m, t - lambdaK[k2, k1, j, m]];

# Second stage mass balance for distribution centers
subject to massBalanceDC2 {j in J, k1 in K_DC, t in 2..n} :
   sum {k2 in K, m in M : (k1, k2) in P2DC_LINKS} F [k1, k2, j, m, t] + sum {r in R, m in M : (k1, r) in DC2R_LINKS} S [k1, r, j, m, t] = I [k1, j, t-1] - I [k1, j, t] + sum {k2 in K, m in M : (k2, k1) in P2DC_LINKS} F [k2, k1, j, m, t - lambdaK[k2, k1, j, m]];

# Second stage mass balance for customers (The connecting constraints First and Second stage)
subject to massBalanceC2 {r in R, j in J, t in T}:
   sum {k in K, m in M : (k, r) in PDC2R_LINKS and t > lambdaR[k, r, j, m]} S [k, r, j, m, t - lambdaR[k, r, j, m]] + SF[r, j, t] = d[r, j, t] ;

# Second stage capacity constraints
subject to capacity2 {k in K_P, j in J, t in 2..n}:
   W[k, j, t] <= Q[k, j];

# Second stage Minimum inventory constraints
subject to minimumInv2 {k in K, j in J, t in 2..n}:
   I[k, j, t] >= Im[k, j, t];
end;
