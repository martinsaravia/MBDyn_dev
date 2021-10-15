%=========================================================================
%
%                  BLADE ELEMENT
%   27-11-2011
% Fuerzas inerciales incrementales, distribuidas totales, internas totales
%==========================================================================
function [X,gdof,m,kt,kg,fie,fde,fme,fae,SE,EE,AV] = scp3llh(FE,X,E,N,D,DM,U,Ui,Vel,Acc,FD,FI,FM,SE,EE,el,AV,RT)

%-----------------------------------------------------------------------------------------------
nxel = 2;      nn1 = E(el,4);         nn2 = E(el,5);   
secid = E(el,3); % Section Properties Location on D matrix
%-----------------------------------------------------------------------------------------------

I0=zeros(3,3);
longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);

% -------------------------   VARIABLES NODALES  ------------------------------------
udofn1=X.CTable(nn1,1:3);       udofn2=X.CTable(nn2,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);       rdofn2=X.CTable(nn2,4:6);  % DOF de Rotaciones
gdof = [ udofn1 rdofn1     udofn2 rdofn2];  

phi1=Ui(rdofn1,X.tstep+1);            phi2=Ui(rdofn2,X.tstep+1);              % Rotación Incremental
x0n1 = N(nn1,2:4)' +  U(udofn1,X.tstep+1);                 x0n2 = N(nn2,2:4)' + U(udofn2,X.tstep+1);                          % Posición actual

idxnn1 = nn1*3-2 : nn1*3;  idxnn2 = nn2*3-2 : nn2*3;

DR1 = RT(idxnn1,1:3,1);       DR2 = RT(idxnn2,1:3,1);   %Incremental Rotation Matrix
T1  =  RT(idxnn1,1:3,2);      T2  = RT(idxnn2,1:3,2);   %Increental Tangential Transformation
 

% Forma vieja de calcular las rotaciones, anda bien. La use para ver problemas en la actualizacion nueva de rotaciones
% phi1=Ui(rdofn1,X.tstep+1);            phi2=Ui(rdofn2,X.tstep+1);              % Rotación Incremental
% DR1 =  expmap(phi1);                 DR2 = expmap(phi2); % Rotacion total  
% T1 =  tangmap (phi1);               T2 = tangmap (phi2); 

en1= DR1 * AV(el*3-2:el*3,1:3,X.tstep);                     en2= DR2 * AV(el*3-2:el*3,4:6,X.tstep); % Current Configuration Triad
AV(el*3-2:el*3,1:6,X.tstep+1)= [ en1  en2 ];  % Guardo las ternas actuales %ASIGNACION DE TERNAS ACTUALES Y MATRICES DE ROTACION ACTUALES

se2n1=skew(en1(:,2));    se3n1=skew(en1(:,3)); %OJO CAMBIO DE SIGNO
se2n2=skew(en2(:,2));    se3n2=skew(en2(:,3)); %OJO CAMBIO DE SIGNO 


% -------------------------   VARIABLES EN ip1, eta=0  ------------------------------------
% Funciones de forma y sus derivadas en ip=0;
N1=0.5;         N2=0.5;           dN1=-1/longel;    dN2=1/longel; 
NI=dia3(N1);         dNI=dia3(dN1);   jac=longel/2; 
     
phiip1=N1*phi1+N2*phi2 ;   Tip1=tangmap(phiip1); 

ei = N1*en1 + N2*en2;        ei2=ei(:,2);       ei3=ei(:,3);
de = dN1*en1 + dN2*en2;      de2=de(:,2);       de3=de(:,3);
dx0 = dN1*x0n1 + dN2*x0n2;  
sei3=skew( N1*en1(:,3)+N2*en2(:,3) ); 
sdx0=skew(dx0);          sde2=skew(de2);  sde3=skew(de3);


% -------------------------------------------------------------------------------------------
%%                                  T A N G E N T      M A T R I C E S 
% --------------------------------------------------------------------------------------------

%% -------------------------------------- TANGENT OF INTERNAL FORCES ---------------------------------------
km=zeros(12,12); kg=zeros(12,12); kt=zeros(12,12); H = zeros (9,18);   G = zeros(18,18);
nips=1; %1 Point Integrations


% Beam Strains and Stresses
EE(1,el,X.tstep+1) = 0.5*( dx0'*dx0 - EE(1,el,1) ); % Axial Strain
EE(2,el,X.tstep+1) = dx0'*de3 - EE(2,el,1);  
EE(3,el,X.tstep+1) = dx0'*de2 - EE(3,el,1); 
EE(4,el,X.tstep+1) = dx0'*ei2 - EE(4,el,1); 
EE(5,el,X.tstep+1) = dx0'*ei3 - EE(5,el,1); 
EE(6,el,X.tstep+1) = de2'*ei3 - EE(6,el,1);
EE(7,el,X.tstep+1) = de2'*de2 - EE(7,el,1);
EE(8,el,X.tstep+1) = de3'*de3 - EE(8,el,1);
EE(9,el,X.tstep+1) = de2'*de3 - EE(9,el,1);

% TENSIONES VIGAS ANTES DE ACTUALIZAR PARA LA GEOMETRICA
SE(1:9,el,X.tstep+1) = D(:,:,secid) * EE(1:9,el,X.tstep+1);  % Total Beam Stresses
Nax1=SE (1,el,X.tstep+1);    M2=SE (2,el,X.tstep+1);     M3=SE (3,el,X.tstep+1);     Q2=SE (4,el,X.tstep+1);  Q3=SE(5,el,X.tstep+1);    
M1=SE(6,el,X.tstep+1);       P2=SE(7,el,X.tstep+1);      P3=SE(8,el,X.tstep+1);      P23=SE(9,el,X.tstep+1);  % Rename Beam Stresses


% H Matrix - Evaluated at integration point ip=0
H(1,1:3) = dx0';
H(2,1:3) = de3';                                                                         H(2,16:18) = dx0';
H(3,1:3) = de2';                                                  H(3,13:15) = dx0';
H(4,1:3) = ei2';    H(4,7:9) = dx0';
H(5,1:3) = ei3';                         H(5,10:12) = dx0';
                                         H(6,10:12) = de2';       H(6,13:15) = ei3';
                                                                  H(7,13:15) = 2*de2';
                                                                                         H(8,16:18) = 2*de3' ;
                                                                  H(9,13:15) = de3' ;    H(9,16:18) = de2' ;
                        
% B Matrix - B = [Bj   Bj]
B = [       dNI                I0                -dNI                I0           ;
            I0                 NI                 I0                 NI           ;    
            I0           N1 * se2n1' * T1'        I0           N2 * se2n2' * T2'    ;
            I0           N1 * se3n1' * T1'        I0           N2 * se3n2' * T2'    ;     
            I0          dN1 * se2n1' * T1'        I0          dN2 * se2n2' * T2'    ;
            I0          dN1 * se3n1' * T1'        I0          dN2 * se3n2' * T2'    ];
 

% A Matrix arising from linealization of Virtual Strains xiT( f_cross( en1(:,3),dx0 ),phi1)
% BE CAREFULL - AAAn MUST BE FED WITH TOTAL BEAM STRESSES
AAAn1 =( M2*dN1+Q3*N1 )*( xiT( f_cross( en1(:,3),dx0 ),phi1)+T1*sdx0*se3n1*T1' ) + ...
        ( M3*dN1+Q2*N1 )*( xiT( f_cross( en1(:,2),dx0 ),phi1)+T1*sdx0*se2n1*T1' ) + ... 
        M1 * ( dN1*     (xiT( f_cross( en1(:,2),ei3 ),phi1)+T1*sei3*se2n1*T1') + N1*( xiT(f_cross(en1(:,3),de2),phi1)+T1*sde2*se3n1*T1') ) + ...
        2*P2 * dN1 * (xiT( f_cross(en1(:,2),de2) ,phi1)+T1*sde2*se2n1*T1') + ... 
        2*P3 * dN1 * (xiT( f_cross(en1(:,3),de3) ,phi1)+T1*sde3*se3n1*T1') + ... 
        P23 * dN1 * ( (xiT( f_cross(en1(:,2),de3 ),phi1)+T1*sde3*se2n1*T1') + (xiT( f_cross(en1(:,3),de2) ,phi1)+T1*sde2*se3n1*T1') ); 
  
AAAn2 = ( M2*dN2+Q3*N2 )*( xiT( f_cross(en2(:,3),dx0) ,phi2)+T2*sdx0*se3n2*T2' ) + ... 
      ( M3*dN2+Q2*N2 )*( xiT( f_cross(en2(:,2),dx0) ,phi2)+T2*sdx0*se2n2*T2' ) + ...
        M1 * ( dN2*(xiT( f_cross(en2(:,2),ei3) ,phi2) + T2*sei3*se2n2*T2') + N2*(xiT( f_cross(en2(:,3),de2) ,phi2)+T2*sde2*se3n2*T2') ) + ...
        2*P2 * dN2 * (xiT( f_cross(en2(:,2),de2) ,phi2)+T2*sde2*se2n2*T2') + ...
        2*P3 * dN2 * (xiT( f_cross(en2(:,3),de3) ,phi2)+T2*sde3*se3n2*T2') + ... 
        P23 * dN2 * ( (xiT( f_cross(en2(:,2),de3) ,phi2)+T2*sde3*se2n2*T2') + (xiT( f_cross(en2(:,3),de2) ,phi2)+T2*sde2*se3n2*T2') ); 
    

    
    
% G Matrix
G(1:3,1:3) = dia3(Nax1);                                G(1:3,7:9)=dia3(Q2);         G(1:3,10:12)=dia3(Q3);       G(1:3,13:15)=dia3(M3);        G(1:3,16:18)=dia3(M2);
                           G(4:6,4:6) = AAAn1+AAAn2;
G(7:9,1:3)=dia3(Q2);
G(10:12,1:3)=dia3(Q3);                                                                                          G(10:12,13:15)=dia3(M1); 
G(13:15,1:3)=dia3(M3);                                                              G(13:15,10:12)=dia3(M1);     G(13:15,13:15)= dia3(2*P2);   G(13:15,16:18)= dia3(P23);
G(16:18,1:3)=dia3(M2);                                                                                          G(16:18,13:15)= dia3(P23);    G(16:18,16:18)= dia3(2*P3);


% Element Tangent Stifness Matrix 

sw = int1d (nips);
for i=1:nips
    if strcmp(FE.step{1},'buckle')==1
        km = sw(2,i) * (B' * H' * D(:,:,secid) * H * B) * jac;
        kg = sw(2,i) * (B' * G * B) * jac;
        kt = km + kg;
    else
        kt = sw(2,i) * (B' * ( H'*D(:,:,secid)*H + G ) * B) * jac;
    end
end


%%  -------------------------------------- TANGENT OF EXTERNAL FORCES --------------------------------------- 
kl=zeros(12,12);   sf=zeros(6,6*nxel);   nips=1;



% if nn2 == 51
%    
%    Dmsn = N2 * FE.sforce(1,5:7)'; %concentrated spatial incremental moment
%    
%    lip1 =  [ I0             I0         ;
%              I0       xiT(Dmsn,phiip1) ];  % Only for 1 point integration
% 
%     sw = int1d (nips); % Numerical Integration
%     for i=1:nips
%         for j=1:6
%             sf(j,j) = ( 1 - sw(1,i) ) / 2;     sf(j,j+6) = ( 1 + sw(1,i) ) / 2;
%         end  
%         kl = sw(2,i) * (sf' * lip1 * sf) * jac   + kl; % Concentrated spatial (fixed) Moment Term
%     end
% 
% 
%     
%     
% %     kl = zeros(12,12); kl(10:12,10:12) = xiT(Dmsn,phi2);
% %     
% %     el
% %     kl
% %     max(max(kl))
% %     
% end

%%  -------------------------------------- TANGENT OF INERTIAL FORCES ---------------------------------------
% NOT CONSIDERED


%%  -------------------------------------- TANGENT KERNEL  ---------------------------------------

kt = kt + kl;



%%  -------------------------------------- TANGENT MASS MATRIX  ---------------------------------------    
m=zeros(12,12);   

if strcmp(FE.step{1,3},'DYNAMIC')==1
    sf=zeros(6,6*nxel);   nips=2;   mip1=zeros(6,6);  % 2 Points Integration
    DMel=DM(:,:,secid);      mip1(1:3,1:3) = DMel(1:3,1:3);
    mip1(4:6,4:6) = Tip1' * DMel(4:6,4:6) * Tip1; % m en punto de integracion eta=0
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
nn0 = E(el,7);
udofn0=X.CTable(nn0,1:3);  rdofn0=X.CTable(nn0,4:6);  

idxnn0 = nn0*3-2 : nn0*3;

DR0 = RT(idxnn0,1:3,1); %

En1  = AV(el*3-2:el*3,1:3,1);          En2=AV(el*3-2:el*3,4:6,1);
    
u0n1 = U(udofn1,X.tstep+1);            u0n2 = U(udofn2,X.tstep+1);          u0n0 = U(udofn0,X.tstep+1);

X0n1 = N(nn1,2:4)' ;                   X0n2 = N(nn2,2:4)' ;                 X0n0 = N(nn0,2:4)'  ;

X0n1r = X0n1-X0n0 ;                    X0n2r = X0n2-X0n0 ;     
x0n1 = X0n1  +  u0n1;                  x0n2 = X0n2  + u0n2;           x0n0 = X0n0 + u0n0;     % Posición actual

% Small Displacement seen by Local Observer
R1 = RT(idxnn1,1:3,3);       R2 = RT(idxnn2,1:3,3);  R0 = RT(idxnn0,1:3,3);  
x1r = x0n1 - x0n0;           x2r = x0n2 - x0n0; 

u0n1r = x1r - R0 * X0n1r;    
u0n2r = x2r - R0 * X0n2r;

DR1r = DR0' *  DR1;                 DR2r = DR0' *  DR2;


    
% X.tn(el*3-2:el*3,1:3,X.tstep+1) = DR1r * X.tn(el*3-2:el*3,1:3,X.tstep);   
% X.tn(el*3-2:el*3,4:6,X.tstep+1) = DR2r * X.tn(el*3-2:el*3,4:6,X.tstep);
% tn1r = X.tn(el*3-2:el*3,1:3,X.tstep+1);                                   tn2r = X.tn(el*3-2:el*3,4:6,X.tstep+1);  

tn1r = ( R0'*R1 ) * En1; 
tn2r = ( R0'*R2 ) * En2;
    
du0r = dN1 * u0n1r + dN2 * u0n2r ;
dX0r = dN1 * X0n1r + dN2 * X0n2r ;

Ei = N1*En1  + N2*En2;         dE = dN1*En1  + dN2*En2;  
Ei1=Ei(:,1);    Ei2=Ei(:,2);   Ei3=Ei(:,3);    dE2=dE(:,2);   dE3=dE(:,3);
%     dX0r = Ei1;

tir = N1 * tn1r + N2 * tn2r;       
dtr = dN1 * tn1r + dN2 * tn2r;
t2ir=tir(:,2);   t3ir=tir(:,3);
dt2r=dtr(:,2);   dt3r=dtr(:,3);
    

Del = D(1:6,1:6,secid);


TEE = zeros (6,1);

TEE(1) = du0r' *  (R0 * dX0r ) ;
TEE(2) = dt3r' *  dX0r  ;  
TEE(3) = dt2r' *  dX0r   ; 
TEE(4) = (R0 * t2ir)' * du0r + t2ir' * dX0r ; 
TEE(5) = (R0 * t3ir)' * du0r + t3ir' * dX0r ; 
TEE(6) = dt2r' *  Ei3  ;
    
if el==1
   EE(1:6,el,X.tstep+1)
   TEE
end

SET = SE(1:9,el,X.tstep+1);
SET(1:6,1) = Del * TEE;  % Total Beam Stresses    


%SEtemp = Del * EE(1:6,el,X.tstep+1);

fie = zeros(12,1); sf=zeros(1,6*nxel); nips = 1;    sw = int1d (nips);

for i=1:nips
%     for j=1:6
%         sf(1,j) = ( 1 - sw(1,i) ) / 2;     sf(1,j+3) = ( 1 + sw(1,i) ) / 2;
%     end  
    fie = sw(2,i) * ( (H*B)' * SET) * jac + fie;     % Total Internal Force
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







    
    
    