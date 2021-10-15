%% ========================================================================
%                   TRIADS FOR BEAM ELEMENT
% NOTA: 27/10/11 MODIFICA RUTINA DE CHANDRU POR RUTINA NUEVA  
%
%        THIRD NODE ALIGNS THE Z AXIS
% in: elemento y nodo
% out: matriz de cosenos directores en COLUMNAS en el centro del elemento

function [A]=triads(N,E,el)

%--------------------------------------------------------------------------
nn1=E(el,4);  nn2=E(el,5);  nn3=E(el,6);  % Numero de nodo del elemento
xn1 = N(nn1,2:4);  xn2 = N(nn2,2:4);  xn3 = N(nn3,2:4);    % Node Positions
%--------------------------------------------------------------------------

d1 = (xn2-xn1)'/norm(xn2-xn1); % Director in direction 1

dt = (xn3-xn1)';  %Temporary director (to define the plane)

d2 = f_cross(dt,d1);  d2=d2/norm(d2);

d3 = f_cross(d1,d2);  d3=d3/norm(d3);   % 

A = [d1 d2 d3];

