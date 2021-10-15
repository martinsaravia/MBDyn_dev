%%==============================================================================================
%
%                          CALCULO DE MATRIZ ANTISIMETRICA
%    IN:  vectores
%    OUT: matrices antisimetricas
%==============================================================================================

function skew=skew(vector)
skew=[     0         -vector(3)    vector(2);
       vector(3)         0        -vector(1);
      -vector(2)     vector(1)        0     ];    