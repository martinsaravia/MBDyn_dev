%=========================================================================
%
%                  BLADE ELEMENT
%   27-11-2011
% Fuerzas inerciales incrementales, distribuidas totales, internas totales
%==========================================================================
function [X,gdof,m,kt,kg,fie,fde,fme,SE,EE,AV] = scb2ll(FE,X,E,N,D,DM,U,Vel,Acc,SE,EE,el,AV)

%-----------------------------------------------------------------------------------------------
nxel = 2;      nn1 = E(el,4);         nn2 = E(el,5);   
secid = E(el,3); % Section Properties Location on D matrix
Del = D(1:6,1:6,secid);
%-----------------------------------------------------------------------------------------------

I0=zeros(3,3);
longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);

% -------------------------   VARIABLES NODALES  ------------------------------------
udofn1=X.CTable(nn1,1:3);       udofn2=X.CTable(nn2,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);       rdofn2=X.CTable(nn2,4:6);  % DOF de Rotaciones
gdof = [ udofn1 rdofn1     udofn2 rdofn2];  

phin1= U(rdofn1,X.tstep+1);            phin2= U(rdofn2,X.tstep+1);              % Rotación Incremental
u0n1 = U(udofn1,X.tstep+1);            u0n2 = U(udofn2,X.tstep+1);
x0n1 = N(nn1,2:4)' +  u0n1;            x0n2 = N(nn2,2:4)' + u0n2;                          % Posición actual
En1  = AV(el*3-2:el*3,1:3,1);          En2=AV(el*3-2:el*3,4:6,1);
AV(el*3-2:el*3,1:6,X.tstep+1) = AV(el*3-2:el*3,1:6,X.tstep);  % Actualizacion de ternas
% -------------------------   VARIABLES EN ip1, eta=0  ------------------------------------
% Funciones de forma y sus derivadas en ip=0;
N1=0.5;         N2=0.5;           dN1=-1/longel;    dN2=1/longel; 
NI=dia3(N1);         dNI=dia3(dN1);   jac=longel/2; 
     
Ei =  N1*En1  + N2*En2;     Ei1=Ei(:,1);   Ei2=Ei(:,2);   Ei3=Ei(:,3);
du0  = dN1*u0n1  + dN2*u0n2;  
dphi = dN1*phin1 + dN2*phin2;
dx0  = dN1*x0n1  + dN2*x0n2; 
phii  = N1*phin1  + N2*phin2; 

%% -------------------------------------- TANGENT OF INTERNAL FORCES ---------------------------------------
km=zeros(12,12); kg=zeros(12,12); kt=zeros(12,12); H = zeros (6,9); 
nips=1; %1 Point Integrations

% Beam Strains and Stresses
EE(1,el,X.tstep+1) = Ei1'  *  du0 ;
EE(2,el,X.tstep+1) = dphi' *  Ei2 ;  
EE(3,el,X.tstep+1) = dphi' * -Ei3 ; 
EE(4,el,X.tstep+1) = Ei2' * du0 - phii' * -Ei3 ; 
EE(5,el,X.tstep+1) = Ei3' * du0 + phii' *  Ei2 ; 
EE(6,el,X.tstep+1) = dphi' *  Ei1 ;

% ESFUERZOS
SE(1:6,el,X.tstep+1) = Del * EE(1:6,el,X.tstep+1);  % Total Beam Stresses
    

% H Matrix - Evaluated at integration point ip=0
H(1,1:3) = Ei1';
                                           H(2,7:9) =  Ei2';
                                           H(3,7:9) = -Ei3';
H(4,1:3) = Ei2';       H(4,4:6) = -Ei3';
H(5,1:3) = Ei3';       H(5,4:6) =  Ei2';
                                            H(6,7:9) =  Ei1';
                   
                                      
% B Matrix - B = [Bj   Bj]
B = [  dNI       I0       -dNI       I0      ;
       I0        NI        I0        NI      ;    
       I0       dNI        I0       -dNI     ];
 
%%  -------------------------------------- STIFFNESS MATRIX  ---------------------------------------
sw = int1d (nips);
for i=1:nips
    kt = sw(2,i) * (B' * ( H'*Del*H ) * B) * jac;
end


%%  -------------------------------------- TANGENT MASS MATRIX  ------------------------------------
m=zeros(12,12);   

if strcmp(FE.step{1,3},'DYNAMIC')==1
    sf=zeros(6,6*nxel);   nips=2;   mip1=zeros(6,6);  % 2 Points Integration
    DMel=DM(:,:,secid);      
    mip1(1:3,1:3) = DMel(1:3,1:3);
    mip1(4:6,4:6) = DMel(4:6,4:6) ; % m en punto de integracion eta=0
    sw = int1d (nips); % Numerical Integration
    for i=1:nips
        for j=1:6
            sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
        end   
        m = sw(2,i) * (sf' * mip1 * sf) * jac   + m;
    end
end

% -------------------------------------------------------------------------------------------
%%                               F O R C E    V E C T O R S 
% --------------------------------------------------------------------------------------------

%% -------------------------------- TOTAL INTERNAL FORCES  -----------------------------------------
fie = zeros(12,1); sf=zeros(1,6*nxel); nips = 1;    sw = int1d (nips);
for i=1:nips
%     for j=1:6
%         sf(1,j) = ( 1 - sw(1,i) ) / 2;     sf(1,j+3) = ( 1 + sw(1,i) ) / 2;
%     end  
    fie = sw(2,i) * ( (H*B)' * SE(1:6,el,X.tstep+1)) * jac + fie;     % Total Internal Force
end


%% --------------------------------  INERTIA FORCES  -----------------------------------------
% Ojo, las fuerzas inerciales son incrementales
fme=zeros(12,1);
if strcmp(FE.step{1,3},'DYNAMIC')==1
    DAccel    = Acc(gdof,X.tstep+1) - Acc(gdof,X.tstep);       
    DVelel = Vel(gdof,X.tstep+1) - Vel(gdof,X.tstep);
    fmacc = m * DAccel ; % Linear part of INCREMENTAL inertia forces with respect to accelerations
    fme =  fmacc + fme;  % Incremental Inertial Force
    
%     fme =  m * Acc(gdof,X.tstep+1);  % TOTAL Inertial Force
end



%% -------------------------------- DISTRIBUTED ELEMENT FORCES  -----------------------------------------
% notar que las fuerzas nodales no se aplican en la rutina del elemento, solo se aplican las fuerzas distribuidas, de cuerpo, etc.
fde = zeros(12,1);   
if isfield(FE,'bforce')==1 % GRAVITY Loads - Assuming that they act in the inertia centroid (no tangent terms because of inertia excentricity)
    sf=zeros(6,6*nxel); fgravel=zeros(6,1);  nips = 1;   sw = int1d (nips);   % Numerical Integration Points
    fgravel(1:3,1) = FE.bforce(1:3)' * X.Dinfo{secid,5} ;  % Element weigth
    for i=1:nips
        for j=1:6
            sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
        end  
        fde = sw(2,i) * sf' * fgravel  * jac + fde;     % Total Internal Force
    end
end










    
    
    