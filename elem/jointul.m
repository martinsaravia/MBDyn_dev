%==========================================================================
%                                JOINTS
%
%  
%    elid 99 = REVOLUTE JOINT -   30/10/2010
%    axis = axialnode1 axialnode2
%==========================================================================

function [m,kt,kg,fie,fde,fme,gdof,X,AV,RT] = jointul (X,E,elid,el,U,Ui,Vel,lmdof,AV,RT)

%-----------------------------------------------------------------------------------------------
nn1 = E(el,4); % SLAVE NODE      
nn2 = E(el,5); % MASTER NODE
rotaxis = E(el,7);
%-----------------------------------------------------------------------------------------------
     
p=1e8;   k=10000;  OO=zeros(3,3);

gdof=[ X.CTable(nn1,1:6) X.CTable(nn2,1:11)];
udofn1=X.CTable(nn1,1:3);          udofn2=X.CTable(nn2,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);          rdofn2=X.CTable(nn2,4:6);  % DOF de Rotacion

% TRIADS UPDATE
idxnn1 = nn1*3-2 : nn1*3;  idxnn2 = nn2*3-2 : nn2*3;
DR1 = RT(idxnn1,1:3,1);       DR2 = RT(idxnn2,1:3,1);  %Incremental Rotation Matrix
T1  =  RT(idxnn1,1:3,2);      T2  = RT(idxnn2,1:3,2);  %Increental Tangential Transformation
enA= DR1 * AV(el*3-2:el*3,1:3,X.tstep);                     enB= DR2 * AV(el*3-2:el*3,4:6,X.tstep); % Current Configuration Triad
AV(el*3-2:el*3,1:6,X.tstep+1) = [ enA  enB ];  % Guardo las ternas actuales %ASIGNACION DE TERNAS ACTUALES Y MATRICES DE ROTACION ACTUALES


% LOOP FOR CONSTRAINTS GRADIENT MATRIX
if elid==90;   % Spherical Joint
    
    dofs = 15;  % dofs + lagrange multipliers
    fde = zeros(dofs,1) ; % No gravity in Joint    
    fme = zeros(dofs,1) ; % No inertia forces in joint
    kg = zeros(dofs,dofs);
    m = zeros (dofs,dofs);
    fde=zeros(dofs,1); %No distributed force
    lms  = Ui ( lmdof , X.tstep+1) ;                     % LAGRANGE MULTIPLIERS VECTOR
    cnst = Ui( udofn1,X.tstep+1) - Ui(udofn2,X.tstep+1);   % DISPLACEMENT CONSTRAINTS VECTOR
    
    I1=eye(3,3);    I0=zeros(3,3);
      
    B = [    I1     I0     -I1      I0  ]; %5x12
      
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(3,3)  ];           

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ];       % CONSTRAINT FORCES
    
end


if elid==91; %REVJUL  % Cylindrical Joint
    
    efic = 0.5;
    cvisc = 0*21.7E6/efic;
    
    
    
    % Global Rotation Axis
    if rotaxis==1; rotplane=[2 3]; end;  if rotaxis==2; rotplane=[1 3]; end ;  if rotaxis==3; rotplane=[1 2]; end    
    
    
    dofs=17;  % dofs + lagrange multipliers
    
    fde = zeros(dofs,1) ; % No gravity in Joint    
    fme = zeros(dofs,1) ; % No inertia forces in joint
    kg = zeros(dofs,dofs);
    m = zeros (dofs,dofs);
    
    lms = Ui ( lmdof , X.tstep+1) ;                         % LAGRANGE MULTIPLIERS VECTOR
    cnst13 = Ui( udofn1,X.tstep+1) - Ui(udofn2,X.tstep+1);  % DISPLACEMENT CONSTRAINTS VECTOR
    cnst4 = enA(:,rotaxis)' * enB(:,rotplane(1));           % ROTATIONAL CONTRAINT 1
    cnst5 = enA(:,rotaxis)' * enB(:,rotplane(2));           % ROTATIONAL CONTRAINT 2
    cnst = [cnst13 ; cnst4 ; cnst5];                        % ALL CONSTRAINTS VECTOR 
    
    I1=eye(3,3);    I0=zeros(3,3); 
    
    axis=enA(:,rotaxis); plane1=enB(:,rotplane(1)); plane2=enB(:,rotplane(2));
    
    DR1A = ( skew(axis)   *  plane1 )' * T1;                 DR1B = ( skew(plane1) *  axis )'   * T2;
    DR2A = ( skew(axis)   *  plane2 )' * T1;                 DR2B = ( skew(plane2) *  axis )'   * T2;  
    
    B = [    I1      I0        -I1      I0   ;
          [0 0 0]    DR1A    [0 0 0]   DR1B  ;
          [0 0 0]    DR2A    [0 0 0]   DR2B ]; %5x12
      
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(5,5)  ];           

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ] ;       % CONSTRAINT FORCES
    
    
%     % VISCOUS FORCES    
    fve = zeros(12,1);     
    Dphii = Ui(rdofn1,X.tstep+1) - Ui(rdofn2,X.tstep+1);
    
% % FORMA VIEJA NO CONSISTENTE
%     Dvel = Vel(rdofn1,X.tstep+1) - Vel(rdofn2,X.tstep+1); % Total Rotational Velocity Difference between nodes    
%     rvel = (Dvel' * axis) * axis; 
%     vforce = [0;0;0; cvisc*rvel];
         

    % SO3 Matrices
    TJ = tangmap (Dphii);  % Joint Incremental Tangent Matrix
    DRJ = expmap (Dphii);  % Joint Incremental Rotation Matrix
    RRJ = expmap (U(rdofn1,X.tstep) - U(rdofn2,X.tstep)); % Joint Reference Rotation Matrix
    RJ = RRJ*DRJ;

    Cv = cvisc * [   T1*TJ            T1*-TJ;
                   -RJ'*T2*TJ        RJ'*T2*TJ]; % forma doble mas lenta Cv2 = [ T1 ; -RJ'*T2 ] * cvisc * [ TJ   -TJ ]
    
    
    fvn = Cv * [ Vel(rdofn1,X.tstep+1) ; Vel(rdofn2,X.tstep+1)];  % Disipative moment vector
    
    fve(4:6,1) = fvn(1:3,1);     % Disipative internal generalized force 1
    fve(10:12,1) = fvn(4:6,1);   % Disipative internal generalized force 2
    
    fie(1:12,1) = fie(1:12,1) + fve;  % Internal Forces Vector
    
    % Info
    OmegaJ = TJ * (Vel(rdofn1,X.tstep+1) - Vel(rdofn2,X.tstep+1));
    
    X.omega = OmegaJ;
    power =efic *  fvn'  * [ Vel(rdofn1,X.tstep+1) ; Vel(rdofn2,X.tstep+1)] ; 
    X.aero(1,8:9,X.tstep+1) = [norm(OmegaJ) power];  % Info table

end


if elid==99;   % revolute joint
          
    axisM=2;  axisS=[1 3] ;
    
    dofs=17;  % dofs + lagrange multipliers
    fde = zeros(dofs,1) ; % No gravity in Joint    
    fme = zeros(dofs,1) ; % No inertia forces in joint
    kg = zeros(dofs,dofs);
    m = zeros (dofs,dofs);
    
    lms = Ui ( lmdof , X.tstep+1) ;                     % LAGRANGE MULTIPLIERS VECTOR
    cnst13 = Ui( udofn1,X.tstep+1) - Ui(udofn2,X.tstep+1);   % DISPLACEMENT CONSTRAINTS VECTOR
    cnst4 = enA(:,axisM)' * enB(:,axisS(1));           % ROTATIONAL CONTRAINT 1
    cnst5 = enA(:,axisM)' * enB(:,axisS(2));           % ROTATIONAL CONTRAINT 2
    cnst = [cnst13 ; cnst4 ; cnst5];               % ALL CONSTRAINTS VECTOR 
    
    I1=[1 0 0;0 1 0;0 0 1];    I0=zeros(3,3);
    
    DR1A = T1' *  skew(enA(:,axisM)) * enB(:,axisS(1))  ;    DR1B = T2' *  skew(enB(:,axisS(1))) * enA(:,axisM) ;   
    
    DR2A = T1' *  skew(enA(:,axisM)) * enB(:,axisS(2))  ;    DR2B = T2' *  skew(enB(:,axisS(2))) * enA(:,axisM) ;  
    
    B = [    I1     I0       -I1      I0  ;
          [0 0 0]  DR1A'   [0 0 0]   DR1B' ;
          [0 0 0]  DR2A'   [0 0 0]   DR2B' ]; %5x12
    
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(5,5)  ];           

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ];       % CONSTRAINT FORCES

end


    
    
    
    