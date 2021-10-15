function [AV]=vmodes(X,M,K,modos)

vect = zeros(X.sdof,modos);

% Autosistema del problema de autovalores standard
A = M \ K; %igual a A=inv(M)*K; 

opts.disp=0;

[vector,lambda]=eigs(A,modos,'SM',opts);      % Resuelvo (entrega todos los modos pero me quedo con algunos)

lambda=real(lambda);           % Elimo la parte imaginaria

[lambda,k]=sort(diag(lambda));  % Ordeno los autovalores en orden ascendente.                                 

% Expresión para frecuencia natural en Hz.
freq=sqrt(lambda(1:modos))/6.283;


vect(X.adof,:)=real(vector(:,k(1:modos)));          % Autovectores para los correspondientes autovales extraidos

% Normalizacion a 1 (MAL, HAY QUE HACERLO PARA CADA COLUMMNA)
vect=vect/(max(max(abs(vect))));

% MANDO fect y freq  A LA MISMA MATRIZ, LA PRIMERA COLUMNA CONTIENE LAS
% FRECUENCIAS Y LAS DEMAS LOS AUTOVECTORES
AV = zeros(length(vect(:,1)), length(freq)+1);
AV(1:modos,1) = freq;
AV(:,2:length(freq)+1) = vect; 