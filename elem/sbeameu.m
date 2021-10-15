%=========================================================================
%
%                               SBEAM ELEMENT MATRICES
%
% Note 1: Total Lagrangian Formulation
% Note 2: Two triads per element (this implies almost always two triads per node, one for each element).
%==========================================================================
function [m,kt,kg,fie,gdof,X] = sbeameu (X,E,N,U,DU,A,el,D,DM)

%-----------------------------------------------------------------------------------------------
nn1 = E(el,4);         nn2 = E(el,5);   
secid = E(el,3); %section id
%-----------------------------------------------------------------------------------------------

I1 = dia(1);
I0 = zeros(3,3);
longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);

% Funciones de forma y sus derivadas en ip=0;
N1=0.5;         N2=0.5;           dN1=-1/longel;    dN2=1/longel;  

% CINEMÁTICA DEL ELEMENTO
udofn1=X.CTable(nn1,1:3);          udofn2=X.CTable(nn2,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);          rdofn2=X.CTable(nn2,4:6);  % DOF de Rotaciones
gdof=[ udofn1 rdofn1               udofn2 rdofn2];


phi1=DU(rdofn1,1);                  phi2=DU(rdofn2,1);
U1=U(udofn1,X.tstep+1);             U2=U(udofn2,X.tstep+1);            % Desplazamiento total
phii = N1*phi1 + N2*phi2;           dphi=dN1*phi1 + dN2*phi2;

X0n1(1:3,1) = N(nn1,2:4)';          X0n2(1:3,1) = N(nn2,2:4)';         % Posicion inicial
x0n1 = X0n1 + U1;                   x0n2 = X0n2 + U2;                      % Posición actual
dx0 = dN1*x0n1 + dN2*x0n2;          dX0 = dN1*X0n1 + dN2*X0n2;    

%        NODO 1                NODO 2  
Ri=expmap(phii);
Ti = tangmap (phii);
Tw = [ I1 I0 I0 I0 ; I0 Ti I0 I0 ; I0 I0 I1 I0; I0 I0 I0 Ti ];
%----------------------------------------------------------
%                        CURVATURA
%----------------------------------------------------------
theta_x = phii(1);     theta_y = phii(2);      theta_z = phii(3);
dtheta_x = dphi(1);    dtheta_y = dphi(2);     dtheta_z = dphi(3);

c = 1 /( 1+(theta_x^2+theta_y^2+theta_z^2)/4 ); 

k=[c*(dtheta_x) + (c*(dtheta_z)*(theta_y))/2 - (c*(dtheta_y)*(theta_z))/2;
  c*(dtheta_y) - (c*(dtheta_z)*(theta_x))/2 + (c*(dtheta_x)*(theta_z))/2;
  c*(dtheta_z) + (c*(dtheta_y)*(theta_x))/2 - (c*(dtheta_x)*(theta_y))/2];

K=skew(k);  % Reemplazo el vector de curvaturas por su antisimétrica.


 %% ===============================================================================================
 %                        A C T U A L I Z A C I O N
 % =================================================================================================

X.e(el*3-2:el*3,1:3)
 
X.e(el*3-2:el*3,1:3) = Ri * X.e(el*3-2:el*3,1:3);

ei1 = X.e(el*3-2:el*3,1);   ei2 = X.e(el*3-2:el*3,2); ei3 = X.e(el*3-2:el*3,3);   
 
X.de(el*3-2:el*3,1) = K * ei1 + Ri * X.de(el*3-2:el*3,1);      % Actualización de la derivada de e  -> de_n+1 = K * e_n+1 + R * de_n
X.de(el*3-2:el*3,2) = K * ei2 + Ri * X.de(el*3-2:el*3,2); 
X.de(el*3-2:el*3,3) = K * ei3 + Ri * X.de(el*3-2:el*3,3);
    
de1=X.de(el*3-2:el*3,1);
de2=X.de(el*3-2:el*3,2);  
de3=X.de(el*3-2:el*3,3);

dA2=X.dA(el*3-2:el*3,2);  
dA3=X.dA(el*3-2:el*3,3);

A2=A(el*3-2:el*3,2);  
A3=A(el*3-2:el*3,3);
  
 %% ===============================================================================================================
 %
 %                                  M A T R I Z   H de vectores base
 %
 % ================================================================================================================
 % Traspongo las matrices x0_int, e y de - 	MEJORAR ESTO
 
 
  b1 = f_cross(ei1,dx0)';       b2 = f_cross(ei2,dx0)';        b3  = f_cross(ei3,dx0)'; % SIGNO POSITIVO (AUNQUE PARA QUE LA MATRIZ DE RIGIDEZ SEA IGUAL A MI ELEMENTO LINEAL DEBO PONERLO NEGATIVO)
  db1 = f_cross(de1,dx0)';      db2 = f_cross(de2,dx0)';       db3 = f_cross(de3,dx0)';
  c2 = (2*f_cross(ei2,de2))';   c3 = (2*f_cross(ei3,de3))';    d = (f_cross(ei2,de3) + f_cross(ei3,de2))';
 
 % Formulacion SARAVIA
 H(1,1:3) = dx0;
 H(2,1:3) = de3;             H(2,4:6) = db3;         H(2,7:9) = b3;    
 H(3,1:3) = de2;             H(3,4:6) = db2;         H(3,7:9) = b2;    
 H(4,1:3) = ei2;             H(4,4:6) = b2;     
 H(5,1:3) = ei3;             H(5,4:6) = b3; 
 
                                                    H(6,7:9) = f_cross(ei2,ei3)'; 

                                                    H(7,7:9) = c2; 
                                                    H(8,7:9) = c3;  
                                                    H(9,7:9) = d; 
 													  												  										  
 %% ===============================================================================================================
 %
 %                                DEFORMACIONES Y TENSIONES INICIALES
 %
 % ================================================================================================================

% % Deformaciones
% ET=[0.5*(dx0'*dx0-1)         dx0'*de3'           dx0'*de2' ... 
%            dx0'*ei2          dx0'*ei3            de2*ei3... 
%            de2*de2'          de3*de3'            de2*de3']';
%        
%        
% Se = D * ET(:,1);        % Esfuerzo Viga
% 
% N1=Se(1);    M2=Se(2);  M3=Se(3);  Q2=Se(4);  Q3=Se(5);  
% T1=Se(6);   P2=Se(7);   P3=Se(8);  P23=Se(9);


EG1 = 0.5*(dx0'*dx0 - dX0'*dX0);
EG2 = dx0'*de3 - dX0'*dA3;
EG3 = dx0'*de2 - dX0'*dA2;
EG4 = dx0'*ei2  - dX0'*A2;
EG5 = dx0'*ei3  - dX0'*A3;
EG6 = de2'*ei3 - dA2'*A3;
EG7 = de2'*de2 - dA2'*dA2;
EG8 = de3'*de3 - dA3'*dA3;
EG9 = de2'*de3 - dA2'*dA3;

EG=[ EG1 EG2 EG3 EG4 EG5 EG6 EG7 EG8 EG9 ]'; % Deformaciones Generalizadas

Se = D(:,:,secid) * EG(:,1);        % Esfuerzo Viga

NN1=Se(1);    M2=Se(2);  M3=Se(3);  Q2=Se(4);  Q3=Se(5);  
T1=Se(6);   P2=Se(7);   P3=Se(8);  P23=Se(9);


%% ================================================================================================
 %
 %                         MATRIZ  G  de SARAVIA
 %
 % =================================================================================================
 
E1 = skew(ei1);        E2 = skew(ei2);     E3 = skew(ei3);
                      dE2 = skew(de2);   dE3 = skew(de3);
dx0 = skew(dx0);
TORS = -T1*E1 + 2*( P2*(dE2*E2+E2*(-dE2)) + P3*(dE3*E3+E3*(-dE3)) ) + P23*( dE2*E3+dE3*E2 + E3*(-dE2)+E2*(-dE3) );
 
G11 =  NN1*I1;     G12 = -(M2*dE3+M3*dE2)-(Q2*E2+Q3*E3);                               G13 = -(M2*E3+M3*E2);
G21 = -G12;       G22 = (M2*dx0*dE3 + M3*dx0*dE2) +(Q2*dx0*E2 + Q3*dx0*E3);           G23 =  M2*dx0*E3 + M3*dx0*E2;
G31 = -G13;       G32 =  G23 + TORS;                                                  G33 =  2*P2*(E2*(-E2)) + 2*P3*(E3*(-E3)) + P23*(E3*(-E2)+E2*(-E3));
 
G = [ G11     	  G12     	 G13;     
     G21     	  G22     	 G23;   
     G31     	  G32        G33];



L(1,1)=dN1;   L(1,7)=dN2;   
L(2,2)=dN1;   L(2,8)=dN2;  
L(3,3)=dN1;   L(3,9)=dN2;  
														  
							L(4,4)=N1;   L(4,10)=N2;   
							L(5,5)=N1;   L(5,11)=N2;  
							L(6,6)=N1;   L(6,12)=N2;														  
															  
														L(7,4)=dN1;   L(7,10)=dN2;   
														L(8,5)=dN1;   L(8,11)=dN2;  
														L(9,6)=dN1;   L(9,12)=dN2;	


%% ------------------------------------------------------------------------
%                   STIFNESS, MASS, AND INTERNAL VECTORS
% ------------------------------------------------------------------------

J=longel/2;      
pint1=0;                                        peso1=2;
pint2=[-0.5773502692 0.5773502692];             peso2=[1 1];

% Integración numérica de la matriz de masa - necesito dos ptos para integrar de forma exacta.
m=zeros(12,12);
NN=zeros(6,12); 

DMel=DM(:,:,secid);
DMel(4:6,4:6) = Ti' * DMel(4:6,4:6) * Ti;

for i=1:2
    eta=pint2(i);
    NN(1,1)=(1-eta)/2;  NN(1,7)=(1+eta)/2;  
    NN(2,2)=(1-eta)/2;  NN(2,8)=(1+eta)/2;
    NN(3,3)=(1-eta)/2;  NN(3,9)=(1+eta)/2;
    NN(4,4)=(1-eta)/2;  NN(4,10)=(1+eta)/2;
    NN(5,5)=(1-eta)/2;  NN(5,11)=(1+eta)/2;
    NN(6,6)=(1-eta)/2;  NN(6,12)=(1+eta)/2;     
    m = peso2(i)*(NN'*DMel*NN)*J + m;
end


% Matriz de Rigidez Material
km = Tw'*(peso1*(L' * H' * D(:,:,secid) * H * L)*J)*Tw;

kg = Tw'*(peso1*(L' * G * L)*J)*Tw;

kt = km + kg;

B=H*L;                  % Ojoooo si pongo directo H y el en fie tengo qe poner L'*H' y no L'*H
fie=peso1*(B'*Se)*J;     % Ojo, tengo que integrar numericamente el esfuerzo viga Se !!!



