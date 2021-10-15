%==========================================================================
%                                JOINTS
%
%  
%    elid 99 = REVOLUTE JOINT -   30/10/2010
%    axis = axialnode1 axialnode2
%==========================================================================

function [m,kt,kg,fie,fbe,X] = joint (X,E,elid,el,U,lmdof)
    
%-----------------------------------------------------------------------------------------------
nn1 = E(el,4);         nn2 = E(el,5);   
%-----------------------------------------------------------------------------------------------

p=1e8;   k=10000;
udofn1=X.CTable(nn1,1:3);          udofn2=X.CTable(nn2,1:3);  % DOF de Desplazamientos
rdofn1=X.CTable(nn1,4:6);          rdofn2=X.CTable(nn2,4:6);  % DOF de Rotaciones             
phi1=U(rdofn1,X.tstep+1);          phi2=U(rdofn2,X.tstep+1);              % node rotation vector
R1 = expmap (phi1);                R2 = expmap (phi2);
T1 = tangmap (phi1);               T2 = tangmap (phi2);

An1 = X.A(el*3-2:el*3,1:3);          An2 = X.A(el*3-2:el*3,4:6);

enA= R1 * An1;                     enB= R2 * An2;

X.e(el*3-2:el*3,1:6)= [ enA  enB ];  % Guardo las ternas actuales

% LOOP FOR CONSTRAINTS GRADIENT MATRIX

if elid==90;   % Spherical Joint
    
    dofs = 15;  % dofs + lagrange multipliers
    lms  = U ( lmdof , X.tstep+1) ;                     % LAGRANGE MULTIPLIERS VECTOR
    cnst = U( udofn1,X.tstep+1) - U(udofn2,X.tstep+1);   % DISPLACEMENT CONSTRAINTS VECTOR
    
    I1=[1 0 0;0 1 0;0 0 1];    I0=zeros(3,3);
      
    B = [    I1     I0     -I1      I0  ]; %5x12
      
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(3,3)  ];           

    kg = zeros(dofs,dofs);

    m = zeros (dofs,dofs);

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ];       % CONSTRAINT FORCES
    fbe = zeros(dofs,1) ; % Massless Joint
end


if elid==91;   % Cylindrical Joint
    
    rotaxis=E(el,5);  
    if rotaxis==1; rotplane=[2 3]; end;  if rotaxis==2; rotplane=[1 3]; end ;  if rotaxis==3; rotplane=[1 2]; end
    
    dofs=17;  % dofs + lagrange multipliers
    lms = U ( lmdof , X.tstep+1) ;                     % LAGRANGE MULTIPLIERS VECTOR
    cnst13 = U( udofn1,X.tstep+1) - U(udofn2,X.tstep+1);   % DISPLACEMENT CONSTRAINTS VECTOR
    cnst4 = enA(:,rotaxis)' * enB(:,rotplane(1));           % ROTATIONAL CONTRAINT 1
    cnst5 = enA(:,rotaxis)' * enB(:,rotplane(2));           % ROTATIONAL CONTRAINT 2
    cnst = [cnst13 ; cnst4 ; cnst5];               % ALL CONSTRAINTS VECTOR 
    
    I1=[1 0 0;0 1 0;0 0 1];    I0=zeros(3,3); 
    
    axis=enA(:,rotaxis); plane1=enB(:,rotplane(1)); plane2=enB(:,rotplane(2));
    
    R1A = ( skew(axis)   *  plane1 )' * T1;                   R1B = ( skew(plane1) *  axis )'   * T2;
    R2A = ( skew(axis)   *  plane2 )' * T1;                   R2B = ( skew(plane2) *  axis )'   * T2;  
    
    B = [    I1      I0       -I1      I0  ;
          [0 0 0]    R1A    [0 0 0]   R1B  ;
          [0 0 0]    R2A    [0 0 0]   R2B ]; %5x12
      
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(5,5)  ];           

    kg = zeros(dofs,dofs);

    m = zeros (dofs,dofs);

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ];       % CONSTRAINT FORCES
    fbe = zeros(dofs,1) ; % Massless Joint
end


if elid==99;   % revolute joint
        
    
%==========================================================================
%                     FORMULACION CON DIRECTORES
%==========================================================================    
    
    axisM=2;  axisS=[1 3] ;
    
    dofs=17;  % dofs + lagrange multipliers
    lms = U ( lmdof , X.tstep+1) ;                     % LAGRANGE MULTIPLIERS VECTOR
    cnst13 = U( udofn1,X.tstep+1) - U(udofn2,X.tstep+1);   % DISPLACEMENT CONSTRAINTS VECTOR
    cnst4 = enA(:,axisM)' * enB(:,axisS(1));           % ROTATIONAL CONTRAINT 1
    cnst5 = enA(:,axisM)' * enB(:,axisS(2));           % ROTATIONAL CONTRAINT 2
    cnst = [cnst13 ; cnst4 ; cnst5];               % ALL CONSTRAINTS VECTOR 
    
    I1=[1 0 0;0 1 0;0 0 1];    I0=zeros(3,3);
    
    R1A = T1' *  skew(enA(:,axisM)) * enB(:,axisS(1))  ;    R1B = T2' *  skew(enB(:,axisS(1))) * enA(:,axisM) ;   
    
    R2A = T1' *  skew(enA(:,axisM)) * enB(:,axisS(2))  ;    R2B = T2' *  skew(enB(:,axisS(2))) * enA(:,axisM) ;  
    
    
    B = [    I1     I0     -I1      I0  ;
      [0 0 0]  R1A'   [0 0 0]   R1B' ;
      [0 0 0]  R2A'   [0 0 0]   R2B' ]; %5x12
      
    
    kt = [ p*(B'*B)      k*B'     ;
             k*B      zeros(5,5)  ];           

    kg = zeros(dofs,dofs);

    m = zeros (dofs,dofs);

    fie = [ B' * (p*cnst+k*lms)   ;   k * cnst ];       % CONSTRAINT FORCES
    fbe = zeros(dofs,1) ; % Massless Joint
end


    
    
    
    