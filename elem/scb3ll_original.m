%=========================================================================
%
%                 LINEAR MULTIBODY BLADE ELEMENT
%   31-10-2012
% Fuerzas inerciales incrementales, distribuidas totales, internas totales
%==========================================================================


function [X,gdof,m,kt,kg,fie,fde,fme,fae,SE,EE,AV] = scb3ll(FE,X,E,N,D,DM,U,Vel,Acc,SE,EE,el,AV,RT)

%-----------------------------------------------------------------------------------------------
nxel = 2;      nn1 = E(el,4);         nn2 = E(el,5);    nn0 = E(el,7);
secid = E(el,3); % Section Properties Location on D matrix
Del = D(1:6,1:6,secid);
%-----------------------------------------------------------------------------------------------

I0=zeros(1,3); I03=zeros(3,3);
longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);


% -------------------------   CTES DE INTERPOLACIÓN EN ip1, eta=0  ------------------------------------
% Funciones de forma y sus derivadas en ip=0;
N1=0.5;         N2=0.5;           dN1=-1/longel;    dN2=1/longel; 
NI=dia3(N1);         dNI=dia3(dN1);   jac=longel/2; 


% -------------------------   VARIABLES NODALES  ------------------------------------
udofn1=X.CTable(nn1,1:3);       udofn2=X.CTable(nn2,1:3);      udofn0=X.CTable(nn0,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);       rdofn2=X.CTable(nn2,4:6);      rdofn0=X.CTable(nn0,4:6);% DOF de Rotaciones

gdof = [ udofn1 rdofn1     udofn2 rdofn2      udofn0 rdofn0];  

% Global Position Vectors
u0n1 = U(udofn1,X.tstep+1);            u0n2 = U(udofn2,X.tstep+1);          u0n0 = U(udofn0,X.tstep+1);
X0n1 = N(nn1,2:4)' ;                   X0n2 = N(nn2,2:4)' ;                 X0n0 = N(nn0,2:4)' ;
x0n1 = N(nn1,2:4)' +  u0n1;            x0n2 = N(nn2,2:4)' + u0n2;           x0n0 = N(nn0,2:4)' + u0n0;     % Posición actual

% Global Rotations
idxnn1 = nn1*3-2 : nn1*3;  idxnn2 = nn2*3-2 : nn2*3;  idxnn0 = nn0*3-2 : nn0*3;
% R1 = RT(idxnn1,1:3,3);         R2 = RT(idxnn2,1:3,3);     R0 = RT(idxnn0,1:3,3); %Total Rotation Matrix
DR1 = RT(idxnn1,1:3,1);        DR2 = RT(idxnn2,1:3,1);    DR0 = RT(idxnn0,1:3,1); %Incremental Rotation Matrix
T1  =  RT(idxnn1,1:3,2);       T2  = RT(idxnn2,1:3,2);     T0  = RT(idxnn0,1:3,2); %Increental Tangential Transformation


% Small Rotation in Global Observer
R1s = R0' *  R1;                    R2s = R0' *  R2;
DR1s = DR0' *  DR1;                 DR2s = DR0' *  DR2;


% Small Displacement seen by Local Observer
u0n1r = R0' * (x0n1 - x0n0) - (X0n1-X0n0);    
u0n2r = R0' * (x0n2 - x0n0) - (X0n2-X0n0);


% TRIAD UPDATE
En1  = AV(el*3-2:el*3,1:3,1);          En2=AV(el*3-2:el*3,4:6,1);

tn1r= R1r * En1;                       tn2r = R2r * En2; 

AV(el*3-2:el*3,1:6,X.tstep+1)= [ tn1r  tn2r ];  

% -------------------------   VARIABLES INTERPOLADAS  ------------------------------------  
du0r = dN1 * u0n1r + dN2 * u0n2r ;
dX0r = dN1 * X0n1r + dN2 * X0n2r ;

Ei =  N1*En1  + N2*En2;     Ei1=Ei(:,1);   Ei2=Ei(:,2);   Ei3=Ei(:,3);

tir = N1 * tn1r + N2 * tn2r;       
dtr = dN1 * tn1r + dN2 * tn2r;

% -------------------------   REDIFINICIONES  ------------------------------------

% t1ir=tir(:,1);   t2ir=tir(:,2);   t3ir=tir(:,3);
% dt1r=dtr(:,1);   dt2r=dtr(:,2);   dt3r=dtr(:,3);
% 
% sE1r = skew(E1r);
% sE2 = skew(En1(:,2));
% sE3 = skew(En1(:,3));
% st2i = skew(tir(:,2));
% st3i = skew(tir(:,3));
% 
% st2n1 = skew(tn1(:,2));
% st2n2 = skew(tn2(:,2));
% st3n1 = skew(tn1(:,3));
% st3n2 = skew(tn2(:,3));
% 
% sdt2 = skew(dt(:,2));
% sdt3 = skew(dt(:,3));
% sdu0r = skew(du0r);

%% -------------------------------------- TANGENT OF INTERNAL FORCES ---------------------------------------
E1r = En1r(:,1);                   
% Beam Strains and Stresses


EE(1,el,X.tstep+1) = E1r'  *  du0r ;
EE(2,el,X.tstep+1) = dt3r' *  E1r ;  
EE(3,el,X.tstep+1) = dt2r' *  E1r ; 
EE(4,el,X.tstep+1) = Ei2' * du0r + t2ir' * E1r ; 
EE(5,el,X.tstep+1) = Ei3' * du0r + t3ir' * E1r ; 
EE(6,el,X.tstep+1) = dt2r' *  Ei3 ;



% ESFUERZOS
SE(1:6,el,X.tstep+1) = Del * EE(1:6,el,X.tstep+1);  % Total Beam Stresses
    

B = [    En1(:,1)' * dN1             I0             En2(:,1)' * dN2             I0        ;
               I0             En1(:,2)' * dN1             I0             En2(:,2)' * dN2  ;
               I0            -En1(:,3)' * dN1             I0            -En2(:,3)' * dN2  ;               
         En1(:,2)' * dN1     -En1(:,3)' * N1        En2(:,2)' * dN2     -En2(:,3)' * N2   ;               
         En1(:,3)' * dN1      En1(:,2)' * N1        En2(:,3)' * dN2      En2(:,2)' * N2   ;
               I0             En1(:,1)' * dN1             I0             En2(:,1)' * dN2  ];

           
Q = [  R0'      I03      I03     I03       -R0'         skew(X0n1r)*T0    ; 
       I03      T0       I03     I03       I03             -T0            ;     
       I03      I03       R0'    I03       -R0'         skew(X0n2r)*T0    ; 
       I03      I03      I03     T0        I03             -T0            ];
%    


%            
% B = [       E1r' * dN1                  I0                   E1r' * dN2                  I0                   du0r' * sE1r * T0     ;
%                I0                -E1r' * st3n1 * dN1              I0                -E1r' * st3n2 * dN2       dt3' * sE1r * T0     ;
%                I0                -E1r' * st2n1 * dN1              I0                -E1r' * st2n2 * dN2       dt2' * sE1r * T0     ;                       
%          En1(:,2)' * dN1         -E1r' * st2n1 * N1          En1(:,2)' * dN2        -E1r' * st2n2 * N2        t2i' * sE1r * T0     ;               
%          En1(:,3)' * dN1         -E1r' * st3n1 * N1          En1(:,3)' * dN2        -E1r' * st3n2 * N2        t3i' * sE1r * T0     ;
%                I0              -En1(:,3)' * st2n1* dN1            I0               -En1(:,3)' * st2n2* dN2            I0             ];   
%            
% 
% Q = [  R0'      I03      I03     I03       -R0'         skew(X0n1r)*T0     ; 
%        I03      T0       I03     I03       I03             -T0           ;     
%        I03      I03       R0'    I03       -R0'         skew(X0n2r)*T0     ; 
%        I03      I03      I03     T0        I03             -T0           ;
%        I03      I03      I03     I03       I03             eye(3)        ];
%        
         
       


%%  -------------------------------------- STIFFNESS MATRIX  ---------------------------------------
kg=zeros(18,18); kt=zeros(18,18);
nips=1; %1 Point Integrations
sw = int1d (nips);
for i=1:nips
    kt = sw(2,i) * (Q' * ( B'*Del*B ) * Q) * jac;
end


%%  -------------------------------------- TANGENT MASS MATRIX  ------------------------------------
m=zeros(18,18);   

if strcmp(FE.step{1,3},'DYNAMIC')==1
    sf=zeros(6,12);   nips=2;   mip1=zeros(6,6);  % 2 Points Integration
    DMel=DM(:,:,secid);      
    mip1(1:3,1:3) = DMel(1:3,1:3);
%     mip1(4:6,4:6) = 0*DMel(4:6,4:6) ; % NO ROTATIONAL MASS
    sw = int1d (nips); % Numerical Integration
    for i=1:nips
        for j=1:6
            sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
        end   
        m(1:12,1:12) = sw(2,i) * (sf' * mip1 * sf) * jac   +  m(1:12,1:12);
    end
end

% -------------------------------------------------------------------------------------------
%%                               F O R C E    V E C T O R S 
% --------------------------------------------------------------------------------------------

%% -------------------------------- TOTAL INTERNAL FORCES  -----------------------------------------
fie = zeros(18,1); sf=zeros(1,6*nxel); nips = 1;    sw = int1d (nips);
for i=1:nips
%     for j=1:6
%         sf(1,j) = ( 1 - sw(1,i) ) / 2;     sf(1,j+3) = ( 1 + sw(1,i) ) / 2;
%     end  
    fie = sw(2,i) * ( (B*Q)' * SE(1:6,el,X.tstep+1)) * jac + fie;     % Total Internal Force
end


%% --------------------------------  INERTIA FORCES  -----------------------------------------
% Ojo, las fuerzas inerciales son incrementales
fme=zeros(18,1);
% if strcmp(FE.step{1,3},'DYNAMIC')==1
%     DAccel    = Acc(gdof,X.tstep+1) - Acc(gdof,X.tstep);       
%     DVelel = Vel(gdof,X.tstep+1) - Vel(gdof,X.tstep);
%     fmacc = m * DAccel ; % Linear part of INCREMENTAL inertia forces with respect to accelerations
%     fme =  fmacc + fme;  % Incremental Inertial Force
%     
% %     fme =  m * Acc(gdof,X.tstep+1);  % TOTAL Inertial Force
% end



%% -------------------------------- DISTRIBUTED ELEMENT FORCES  -----------------------------------------
% notar que las fuerzas nodales no se aplican en la rutina del elemento, solo se aplican las fuerzas distribuidas, de cuerpo, etc.
fde = zeros(18,1);   
if isfield(FE,'bforce')==1 % GRAVITY Loads - Assuming that they act in the inertia centroid (no tangent terms because of inertia excentricity)
    sf=zeros(6,6*nxel); fgravel=zeros(6,1);  nips = 1;   sw = int1d (nips);   % Numerical Integration Points
    fgravel(1:3,1) = FE.bforce(1:3)' * X.Dinfo{secid,5} ;  % Element weigth
    for i=1:nips
        for j=1:6
            sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
        end  
        fde(1:12) = sw(2,i) * sf' * fgravel  * jac + fde(1:12);     % Total Internal Force
    end
end


%% -------------------------------- AERODYNAMIC ELEMENT FORCES  -----------------------------------------
% Aerodynamic forces per unit length
fae = zeros(12,1);   

if X.flag(10)==1 % BEM Aerodynamic forces flag
    sf=zeros(6,6*nxel);  nips = 1;   sw = int1d (nips);   % Numerical Integration Points
     [X,fa] = aerobem (X,FE,E,N,AV,el,U,Vel);
    for i=1:nips
        for j=1:6
            sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
        end  
        fae = sw(2,i) * sf' * fa  * jac + fae;     % Total Internal Force
    end
end


% if el ==2 
%     fie(1,13:18)=0;
%     kt(13:18,13:18)=0;
% end


% kt * U(gdof,X.tstep+1)
% squeeze(EE(1:6,el,X.tstep+1))

end

% 
% Ei2' * du0r + ti2' * R0'*Ei1  
% dx0ir'*Ei2 + du0r'*ti2
% 
% dx0ir'*ti2

    