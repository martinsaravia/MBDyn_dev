%=========================================================================
%
%                 LINEAR MULTIBODY BLADE ELEMENT
%   31-10-2012
% Fuerzas inerciales incrementales, distribuidas totales, internas totales
%==========================================================================


function [X,gdof,m,kt,kg,fie,fde,fme,fae,SE,EE,AV] = scb3ll2(FE,X,E,N,D,DM,U,Ui,Vel,Acc,SE,EE,el,AV,RT)

formu = 'INCREMENTAL';
%-----------------------------------------------------------------------------------------------
nxel = 2;      nn1 = E(el,4);         nn2 = E(el,5);    nn0 = E(el,7);
secid = E(el,3); % Section Properties Location on D matrix
Del = D(1:6,1:6,secid);
%-----------------------------------------------------------------------------------------------

I0=zeros(1,3); I03=zeros(3,3);
longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);


% -------------------------   CTES DE INTERPOLACIÓN EN ip1, eta=0  ------------------------------------
% Funciones de forma y sus derivadas en ip=0;
N1=0.5;              N2=0.5;           dN1=-1/longel;    dN2=1/longel; 
NI=dia3(N1);         dNI=dia3(dN1);    jac=longel/2; 


% -------------------------   VARIABLES NODALES  ------------------------------------
udofn1=X.CTable(nn1,1:3);       udofn2=X.CTable(nn2,1:3);      udofn0=X.CTable(nn0,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);       rdofn2=X.CTable(nn2,4:6);      rdofn0=X.CTable(nn0,4:6);% DOF de Rotaciones

gdof = [ udofn1 rdofn1     udofn2 rdofn2      udofn0 rdofn0];  

% Global Rotations
idxnn1 = nn1*3-2 : nn1*3;  idxnn2 = nn2*3-2 : nn2*3;  idxnn0 = nn0*3-2 : nn0*3;


R1 = RT(idxnn1,1:3,3);          R2 = RT(idxnn2,1:3,3);     R0 = RT(idxnn0,1:3,3); %Total Rotation Matrix
DR1 = RT(idxnn1,1:3,1);         DR2 = RT(idxnn2,1:3,1);    DR0 = RT(idxnn0,1:3,1); %Incremental Rotation Matrix
T0  = RT(idxnn0,1:3,2); %Increental Tangential Transformation

% R1 = RT(idxnn1,1:3,3);                      R2 = RT(idxnn2,1:3,3);                   R0 = RT(idxnn0,1:3,3); %Total Rotation Matrix
% DR1 = expmap(Ui(rdofn1,X.tstep+1));         DR2 = expmap(Ui(rdofn2,X.tstep+1));      DR0 = expmap(Ui(rdofn0,X.tstep+1)); %Incremental Rotation Matrix
% T1  = tangmap(Ui(rdofn1,X.tstep+1));        T2  = tangmap(Ui(rdofn2,X.tstep+1));     T0  = tangmap(Ui(rdofn0,X.tstep+1)); %Increental Tangential Transformation


% Global Position Vectors
if strcmp(formu, 'TOTAL')==1;
    
    T0 = tangmap(U(rdofn0,X.tstep+1));
    
    En1  = AV(el*3-2:el*3,1:3,1);                  En2=AV(el*3-2:el*3,4:6,1);
    
    u0n1 = U(udofn1,X.tstep+1);            u0n2 = U(udofn2,X.tstep+1);          u0n0 = U(udofn0,X.tstep+1);
    X0n1 = N(nn1,2:4)' ;                   X0n2 = N(nn2,2:4)' ;                 X0n0 = N(nn0,2:4)' ;
    X0n1r = X0n1-X0n0 ;                    X0n2r = X0n2-X0n0 ;                 
    x0n1 = N(nn1,2:4)' +  u0n1;            x0n2 = N(nn2,2:4)' + u0n2;           x0n0 = N(nn0,2:4)' + u0n0;     % Posición actual
    
        % Small Displacement seen by Local Observer
    u0n1r = R0' * (x0n1 - x0n0) - (X0n1-X0n0);    
    u0n2r = R0' * (x0n2 - x0n0) - (X0n2-X0n0);
    
%     TL FORMULATION
    R1r = R0' *  R1;                    R2r = R0' *  R2;
    en1= R1 * AV(el*3-2:el*3,1:3,1);                     en2= R2 * AV(el*3-2:el*3,4:6,1); % Current Configuration Triad
    tn1r = R1r * En1;                              tn2r = R2r * En2; 
    
    AV(el*3-2:el*3,1:6,X.tstep+1)= [ en1  en2 ];
    
    du0r = dN1 * u0n1r + dN2 * u0n2r ;
    dX0r = dN1 * X0n1r + dN2 * X0n2r ;

    Ei =  N1*En1  + N2*En2;     Ei1=Ei(:,1);   Ei2=Ei(:,2);   Ei3=Ei(:,3);

    tir = N1 * tn1r + N2 * tn2r;       
    dtr = dN1 * tn1r + dN2 * tn2r;
    t1ir=tir(:,1);   t2ir=tir(:,2);   t3ir=tir(:,3);
    dt1r=dtr(:,1);   dt2r=dtr(:,2);   dt3r=dtr(:,3);
    
    
    EE(1,el,X.tstep+1) = dX0r' *  du0r ;
    EE(2,el,X.tstep+1) = dt3r' *  dX0r ;  
    EE(3,el,X.tstep+1) = dt2r' *  dX0r ; 
    EE(4,el,X.tstep+1) = Ei2' * du0r + t2ir' * dX0r ; 
    EE(5,el,X.tstep+1) = Ei3' * du0r + t3ir' * dX0r ; 
    EE(6,el,X.tstep+1) = dt2r' *  Ei3;

    Q = [  R0'      I03       I03     I03       -R0'           skew(R0'*(x0n1-x0n0))*T0    ; 
           I03      T0        I03     I03       I03                     -T0              ; 

           I03      I03       R0'     I03       -R0'           skew(R0'*(x0n2-x0n0))*T0   ; 
           I03      I03       I03     T0        I03                      -T0             ];
end




if strcmp(formu, 'INCREMENTAL')==1;
    
    
    En1  = AV(el*3-2:el*3,1:3,1);                  En2=AV(el*3-2:el*3,4:6,1);
    
    u0n1 = U(udofn1,X.tstep+1);            u0n2 = U(udofn2,X.tstep+1);          u0n0 = U(udofn0,X.tstep+1);
    X0n1 = N(nn1,2:4)' ;                   X0n2 = N(nn2,2:4)' ;                 X0n0 = N(nn0,2:4)'  ;
    
    X0n1r = X0n1-X0n0 ;                    X0n2r = X0n2-X0n0 ;     
    x0n1 = X0n1  +  u0n1;                  x0n2 = X0n2  + u0n2;           x0n0 = X0n0 + u0n0;     % Posición actual
    
    % Small Displacement seen by Local Observer
      R0ref = RT(idxnn0,1:3,4);  
      R0 = R0ref * DR0;
      
%       EE0 = AV(52*3-2:52*3,4:6,1); 
%       ee0 = expmap(Ui(rdofn0,X.tstep+1)) * AV(52*3-2:52*3,4:6,X.tstep);
%       ee0ref = AV(52*3-2:52*3,4:6,X.tstep);
%       
%       R0ref = ee0ref * EE0'; 
%       R0 = ee0 * EE0';
%       
%       TEMP = R0 - R0ref * DR0
      
      u0n1r = R0' * (x0n1 - x0n0) - (X0n1-X0n0);    
      u0n2r = R0' * (x0n2 - x0n0) - (X0n2-X0n0);

  
      DR1r = DR0' *  DR1;                 DR2r = DR0' *  DR2;
    
%     DR1r = expmap (Ui(rdofn1,X.tstep+1)-Ui(rdofn0,X.tstep+1));  
%     DR2r = expmap (Ui(rdofn2,X.tstep+1)-Ui(rdofn0,X.tstep+1));
    
    
    X.tn(el*3-2:el*3,1:3,X.tstep+1) = DR1r * X.tn(el*3-2:el*3,1:3,X.tstep);   
    X.tn(el*3-2:el*3,4:6,X.tstep+1) = DR2r * X.tn(el*3-2:el*3,4:6,X.tstep);
    
    tn1r = X.tn(el*3-2:el*3,1:3,X.tstep+1);                                   tn2r = X.tn(el*3-2:el*3,4:6,X.tstep+1); 
    
%     en1= DR1 * AV(el*3-2:el*3,1:3,X.tstep);                                   en2= DR2 * AV(el*3-2:el*3,4:6,X.tstep); % Current Configuration Triad 
    
    % PRUEBA
%     en1= DR1 * AV(el*3-2:el*3,1:3,X.tstep);                     en2= DR2 * AV(el*3-2:el*3,4:6,X.tstep);
    en1 = R0 * tn1r; en2 = R0 * tn2r;    
    
    en1 = R0 * En1; en2 = R0 * En2; 
    
    AV(el*3-2:el*3,1:6,X.tstep+1)= [ en1  en2 ];
    
    
%     
      R1r = R0' *  R1;                    R2r = R0' *  R2;
%     en1= R0 * AV(el*3-2:el*3,1:3,1);            en2= R0 * AV(el*3-2:el*3,4:6,1); % Current Configuration Triad
%      tn1r = R1r * En1;                              tn2r = R2r * En2;    
%           AV(el*3-2:el*3,1:6,X.tstep+1)= [ en1  en2 ];
%     
    
    du0r = dN1 * u0n1r + dN2 * u0n2r ;
    dX0r = dN1 * X0n1r + dN2 * X0n2r ;

    
    Ei = N1*En1  + N2*En2;         dE = dN1*En1  + dN2*En2;  
    Ei1=Ei(:,1);    Ei2=Ei(:,2);   Ei3=Ei(:,3);    dE2=dE(:,2);   dE3=dE(:,3);
%     dX0r = Ei1;
    
    tir = N1 * tn1r + N2 * tn2r;       
    dtr = dN1 * tn1r + dN2 * tn2r;
    t2ir=tir(:,2);   t3ir=tir(:,3);
    dt2r=dtr(:,2);   dt3r=dtr(:,3);
    
% BORRAR
% phir1 = Ui(rdofn1,X.tstep+1) - Ui(rdofn0,X.tstep+1); phir2 = Ui(rdofn2,X.tstep+1) - Ui(rdofn0,X.tstep+1);
% dphir = dN1 * phir1  + dN2 * phir2;
% 
% dt3r = skew(dphir) * Ei3; 
% dt2r = skew(dphir) * Ei2; 

 
EE(1,el,X.tstep+1) = dX0r' *  du0r  ;
EE(2,el,X.tstep+1) = dt3r' *  dX0r  ;  
EE(3,el,X.tstep+1) = dt2r' *  dX0r   ; 
EE(4,el,X.tstep+1) = Ei2' * du0r + t2ir' * dX0r ; 
EE(5,el,X.tstep+1) = Ei3' * du0r + t3ir' * dX0r ; 
EE(6,el,X.tstep+1) = dt2r' *  Ei3   ;



Q = [ R0'          I03        I03             I03         -R0'         skew(R0'*(x0n1-x0n0))*T0    ; 
      I03          T0         I03             I03          I03                  -T0                ; 

      I03          I03         R0'            I03         -R0'         skew(R0'*(x0n2-x0n0))*T0    ; 
      I03          I03        I03             T0           I03                  -T0                ]; 
  

    
end


% ESFUERZOS
SE(1:6,el,X.tstep+1) = Del * EE(1:6,el,X.tstep+1);  % Total Beam Stresses
    
sE2 = skew (En1(:,2));  sE3 = skew (En1(:,3)); 


B = [    dX0r' * dN1             I0                      dX0r' * dN2                 I0           ;
    
               I0             dX0r' * sE3' * dN1             I0             dX0r' * sE3' * dN2    ;
               I0             dX0r' * sE2' * dN1             I0             dX0r' * sE2' * dN2    ; 
               
         En1(:,2)' * dN1     dX0r' * sE2' * N1        En2(:,2)' * dN2       dX0r' * sE2' * N2     ;
         En1(:,3)' * dN1     dX0r' * sE3' * N1        En2(:,3)' * dN2       dX0r' * sE3' * N2     ;
         
               I0          En1(:,3)' * sE2' * dN1             I0           En1(:,3)' * sE2' * dN2 ];

   


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
fae = zeros(18,1);   
  





    