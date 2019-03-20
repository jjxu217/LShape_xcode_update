* ENCODING=ISO-8859-1
NAME          CEP.lp
ROWS
 N  FOBJ    
 L  HRSM1   
 L  HRSM2   
 L  HRSM3   
 L  HRSM4   
 L  MAI     
 G  CAPM1   
 G  CAPM2   
 G  CAPM3   
 G  CAPM4   
 G  DEMP1   
 G  DEMP2   
 G  DEMP3   
COLUMNS
    xM1       FOBJ                          2.5
    xM1       HRSM1                          -1
    xM2       FOBJ                         3.75
    xM2       HRSM2                          -1
    xM3       FOBJ                            5
    xM3       HRSM3                          -1
    xM4       FOBJ                            3
    xM4       HRSM4                          -1
    zM1       HRSM1                           1
    zM1       MAI                          0.08
    zM1       CAPM1                           1
    zM2       HRSM2                           1
    zM2       MAI                          0.04
    zM2       CAPM2                           1
    zM3       HRSM3                           1
    zM3       MAI                          0.03
    zM3       CAPM3                           1
    zM4       HRSM4                           1
    zM4       MAI                          0.01
    zM4       CAPM4                           1
    yP1M1     FOBJ                          2.6
    yP1M1     CAPM1                          -1
    yP1M1     DEMP1                         0.6
    yP1M2     FOBJ                          3.4
    yP1M2     CAPM2                          -1
    yP1M2     DEMP1                         0.6
    yP1M3     FOBJ                          3.4
    yP1M3     CAPM3                          -1
    yP1M3     DEMP1                         0.9
    yP1M4     FOBJ                          2.5
    yP1M4     CAPM4                          -1
    yP1M4     DEMP1                         0.8
    yP2M1     FOBJ                          1.5
    yP2M1     CAPM1                          -1
    yP2M1     DEMP2                         0.1
    yP2M2     FOBJ                          2.4
    yP2M2     CAPM2                          -1
    yP2M2     DEMP2                         0.9
    yP2M3     FOBJ                            2
    yP2M3     CAPM3                          -1
    yP2M3     DEMP2                         0.6
    yP2M4     FOBJ                          3.6
    yP2M4     CAPM4                          -1
    yP2M4     DEMP2                         0.8
    yP3M1     FOBJ                            4
    yP3M1     CAPM1                          -1
    yP3M1     DEMP3                        0.05
    yP3M2     FOBJ                          3.8
    yP3M2     CAPM2                          -1
    yP3M2     DEMP3                         0.2
    yP3M3     FOBJ                          3.5
    yP3M3     CAPM3                          -1
    yP3M3     DEMP3                         0.5
    yP3M4     FOBJ                          3.2
    yP3M4     CAPM4                          -1
    yP3M4     DEMP3                         0.8
    sP1       FOBJ                          400
    sP1       DEMP1                           1
    sP2       FOBJ                          400
    sP2       DEMP2                           1
    sP3       FOBJ                          400
    sP3       DEMP3                           1
RHS
    RHS       HRSM1                         500
    RHS       HRSM2                         500
    RHS       HRSM3                         500
    RHS       HRSM4                         500
    RHS       MAI                           100
    RHS       DEMP1                        1500
    RHS       DEMP2                        1500
    RHS       DEMP3                        1500
BOUNDS
 UP bnd       zM1                          2000
 UP bnd       zM2                          2000
 UP bnd       zM3                          3000
 UP bnd       zM4                          3000
ENDATA
