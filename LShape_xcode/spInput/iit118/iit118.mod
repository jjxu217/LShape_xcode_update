#### AMPL Model file for two-stage economic dispatch problem with renew energy resources and deterministic demand
#### The model can be written for a single time period or in a rolling horizon setting depending on the flag ROLLING

param ROLLING;

### Decision epochs
param N;               # Number of time periods in the horizon
param F;               # Number of decision epochs in an hour (frequency)
param hour;            # Hour for which the model is written
set TIME := 1..N;

### Power network elements
## Buses
set BUS;
param busName{BUS}, symbolic; param busRegion{BUS}; param vMin{BUS}; param vMax{BUS};
table bus_data IN "CSV" "../../energySystems/iit118/busData.csv" : 
BUS <- [BUS], busName, busRegion, vMax, vMin;

param busV{BUS} := 100;

## Lines
set LINES;
param lineName{LINES}, symbolic; param fromBus{LINES}; param toBus{LINES}; param circuitID{LINES}; param resistance{LINES}; param reactance{LINES};
param susceptance{LINES}; param flowMax{LINES}; param flowMin{LINES};
table lines_data IN "CSV" "../../energySystems/iit118/lineData.csv" : 
LINES <- [LINES], lineName, fromBus, toBus, circuitID, resistance, reactance, susceptance, flowMax, flowMin;

## Conventional generator
set GEN;
param genName{GEN}, symbolic; param genBus{GEN}; param genMax{GEN}; param genMin{GEN};
param genH0{GEN}; param genH1{GEN}; param genH2{GEN}; param genH3{GEN}; param genH4{GEN}; param genH5{GEN};
param genRampUp{GEN}; param genRampDown{GEN}; param genTimeUp{GEN}; param genTimeDown{GEN};
param genStartCost{GEN}; param genOpCost{GEN}; param prodCost;
table gen_data IN "CSV" "../../energySystems/iit118/genData.csv" : 
GEN <- [GEN], genName ~ Name, genBus ~ bus, genMax ~ capMax, genMin ~ capMin, genH0 ~ heatRate_base, genH1 ~ heatRate_1, genH2 ~ heatRate_2, genH3 ~ heatRate_3, genH4 ~ heatRate_4, genH5 ~ heatRate_5,
genRampUp ~ rampUp, genRampDown ~ rampDown, genTimeUp ~ timeUp, genTimeDown ~ timeDown, genStartCost ~ startCost, genOpCost ~ opCost;

## Renewable generator
set RENEW;
param renewName{RENEW}, symbolic; param renewBus{RENEW}; param renewMax{RENEW}; param renewMin{RENEW};
param renewH0{RENEW}; param renewH1{RENEW}; param renewH2{RENEW}; param renewH3{RENEW}; param renewH4{RENEW}; param renewH5{RENEW};
param renewRampUp{RENEW}; param renewRampDown{RENEW}; param renewTimeUp{RENEW}; param renewTimeDown{RENEW};
param renewStartCost{RENEW}; param renewOpCost{RENEW};
table renew_data IN "CSV" "../../energySystems/iit118/renewData.csv" : 
RENEW <- [RENEW], renewName ~ Name, renewBus ~ bus, renewMax ~ capMax, renewMin ~ capMin, renewH0 ~ heatRate_base, renewH1 ~ heatRate_1, renewH2 ~ heatRate_2, renewH3 ~ heatRate_3, renewH4 ~ heatRate_4, renewH5 ~ heatRate_5, renewRampUp ~ rampUp, renewRampDown ~ rampDown, renewTimeUp ~ timeUp, renewTimeDown ~ timeDown, renewStartCost ~ startCost, renewOpCost ~ opCost;

## Load
set LOAD;
param loadBus{LOAD}; param loadRegion{LOAD}; param Pd{LOAD}; param Qd{LOAD};
table load_data IN "CSV" "../../energySystems/iit118/loadData.csv" : 
LOAD <- [LOAD], loadBus, loadRegion, Pd, Qd;

## Fixed parameters: initial generation, deterministic demand, and renew generation
param gInit{g in GEN} := (genMin[g]+genMax[g])/2;
     
set D, dimen 2;
param Period; param Region; param  realLoad{D}; 
table realLoad_data IN "CSV" "../../energySystems/iit118/realLoad_glpk.csv" :
D <- [Period,Region], realLoad;

set W, dimen 2;
param renewPower{W};
table renewPower_data IN "CSV" "../../energySystems/iit118/renewPower_glpk.csv" :
W <- [Period, RENEW], renewPower;

### Optimization model
## Decision variables
# First-stage 
var g0{g in GEN} <= genMax[g], >= genMin[g];

# Second-stage 
var g1{n in TIME, g in GEN: n <> 1} <= genMax[g], >= genMin[g];
var gs1{n in TIME, g in GEN} >= 0;
var p1{n in TIME, l in LINES} <= flowMax[l], >= flowMin[l];
var th1{n in TIME, b in BUS} <= 3.14/4 >= -3.14/4;
# var th1{n in TIME, b in BUS};
var w1{n in TIME, r in RENEW} >= 0;
	 
## Constraints
# First-stage: Ramp rate limits
RU{g in GEN}: g0[g] - gInit[g] <= genRampUp[g]/F;
RD{g in GEN}: g0[g] - gInit[g] >= -genRampDown[g]/F;
		    
# Second-stage: Ramp rate limits
RU0{g in GEN}: g1[2,g] - g0[g] <= genRampUp[g]/F;
RD0{g in GEN}: g1[2,g] - g0[g] >= -genRampDown[g]/F;

RU1{n in TIME, g in GEN: n > 1 && n <> N}: g1[n+1,g] - g1[n,g] <= genRampUp[g]/F;
RD1{n in TIME, g in GEN: n > 1 && n <> N}: g1[n+1,g] - g1[n,g] >= -genRampDown[g]/F;

# Second-stage: flow balance
PF0{n in TIME, i in BUS: n == 1}: (sum{l in LINES: i == toBus[l]} p1[n,l]) 
	   - (sum{l in LINES: i == fromBus[l]} p1[n,l])
           + (sum{g in GEN: i == genBus[g]} g0[g])
           - (sum{g in GEN: i == genBus[g]} gs1[n,g])
	   + (sum{r in RENEW: i == renewBus[r]} w1[n,r])	  
	   = (sum{d in LOAD: i == loadBus[d]} realLoad[floor(hour + (n-1)/F), busRegion[i]]*Pd[d]);

PF1{n in TIME, i in BUS: n <> 1}: (sum{l in LINES: i == toBus[l]} p1[n,l]) 
	   - (sum{l in LINES: i == fromBus[l]} p1[n,l])
           + (sum{g in GEN: i == genBus[g]} g1[n,g])
           - (sum{g in GEN: i == genBus[g]} gs1[n,g])
	   + (sum{r in RENEW: i == renewBus[r]} w1[n,r])  
	   = (sum{d in LOAD: i == loadBus[d]} realLoad[floor(hour + (n-1)/F), busRegion[i]]*Pd[d]);

FE1{n in TIME, l in LINES}: p1[n,l] - (busV[fromBus[l]]*busV[toBus[l]]/reactance[l])*(th1[n,fromBus[l]] - th1[n,toBus[l]]) = 0;

S0{n in TIME, g in GEN: n == 1}: gs1[n,g] <= g0[g];
S1{n in TIME, g in GEN: n <> 1}: gs1[n,g] <= g1[n,g];

R1{n in TIME, r in RENEW}: w1[n,r] <= renewPower[floor(hour + (n-1)/F), r];

## Objective function
minimize totalCost: sum{g in GEN} 1*g0[g] + sum{n in TIME: n <> 1} (sum{g in GEN} 1*g1[n,g]);


data;
param N := 6;
param F := 6;
param hour := 1;
end;
