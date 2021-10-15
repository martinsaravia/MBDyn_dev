%% NOTES: only one wind card, only one rotor card


function [X,fa] = aerobem (X,FE,E,N,AV,el,U,Vel)

metodo = 'NOLINEAL';

if isempty(FE.wind{1,4}) == 0  % if file reading is active
    wtb = FE.wind {2,1};
    time=X.T(2,X.tstep);
    if length(wtb(:,1)) >=2
        wdata = interp1q( wtb(:,1) , wtb , time);% Interpolate de column
        if isnan(wdata(1))==1;
            wdata = interp1( wtb(:,1) , wtb , time, 'nearest','extrap' ); % For out of table data extrapolate to nearest
        end
    else
        wdata = wtb;       % If number of rows is 1 does not interpolate
        if wdata(1) > time  % Exit if time has not reach the starting time
            fa = zeros(6,1);
            return
        end
    end

end

%-----------------------------------------------------------------------------------------------
setname = FE.rotor{1,3}; % only one rotor allowed (sino tengo que armar algo para acceder al nidx)  
% setsec = strmatch(setname,FE.sec(1:2:X.secs,3),'exact')
nidx = E(el,8)*2-1; % LOCATION OF SECTION FILE

ftb = FE.sec{ nidx,5};     btb = FE.sec{nidx,6};  
 
hubnode = FE.rotor{2,1}(1); % Hub Node Number
joint = FE.rotor{2,1}(2)  ; % Joint element Number

% Axis of rotation (joint) and dofs to measure angular speed
jaxis=E(joint,6);  
mdof= X.CTable(E(joint,4),4:6);   
sdof=X.CTable(E(joint,5),4:6);  

rho = 1.2;   % Air density O

%-----------------------------------------------------------------------------------------------
nn1 = E(el,4);                                  nn2 = E(el,5);                               
udofn1 = X.CTable(nn1,1:3);                     udofn2=X.CTable(nn2,1:3);                  udofhub=X.CTable(hubnode,1:3); % DOF de Desplazamientos
xn1 = N(nn1,2:4)' +  U(udofn1,X.tstep+1);

xn2 = N(nn2,2:4)' + U(udofn2,X.tstep+1);   
xhub = N(hubnode,2:4)' + U(udofhub,X.tstep+1);

xc = 0.5*xn1+0.5*xn2;         xcr = xc - xhub;             % Element Centroid Global and relative  Positions
ri = norm(xn1-xhub);          re = norm(xn2-xhub);     rm=(ri+re)/2;    dr = (re - ri);    % Internal Radius    % External Radius     % Delta Radius

%ELEMENT TRIADS
if strcmp(metodo,'NOLINEAL')
    en1 = AV(el*3-2:el*3,1:3,X.tstep+1);        en2 = AV(el*3-2:el*3,4:6,X.tstep+1); % Current Configuration Triad
    
elseif strcmp(metodo,'LINEAL')
    en1 = AV(el*3-2:el*3,1:3,1);        en2 = AV(el*3-2:el*3,4:6,1); % Current Configuration Triad
end
ei = 0.5 * (en1 + en2); % Current Triad at the element center.

% UNIT NORMAL AND TANGENTIAL VECTORS
en = FE.rotor{2,1}(3:5) ;   en=en'/norm(en);  % Rotor Unit Normal Vector 
et = f_cross(en,xc);        et=et/norm(et);   % TANGENTIAL UNIT VECTOR

% GLOBAL VECTORS
Ea = FE.rotor{2,1}(3:5); 
Ev = FE.rotor{2,1}(6:8);  % Vertical Direction
Et = f_cross(Ea,Ev)';      % Transversal direction



%-----------------------------------------------------------------------------------------------
%                                 INPUT VELOCITES                        
%-----------------------------------------------------------------------------------------------

Wr = Vel(mdof,X.tstep+1) - Vel(sdof,X.tstep+1) ;  ome=norm(Wr) ;  % Angular Speed
Vr = - f_cross (Wr,xcr); % Velocity of the WIND caused by the rotation of the blade, LOOK at the minus sign
VB = 0.5 * ( Vel(udofn1,X.tstep+1) + Vel(udofn2,X.tstep+1) ); % Blade Velocity ojooOOOOOOOOOOOOO CON EL PASO DE TIEMPO ES TSTEP+1
% --------------- WIND VELOCITY ---------------
if isempty(FE.wind{1,4}) == 0      
   VW = iecspeed (FE,X,wdata,xcr,xhub,Ea,Ev,Et);            % Obtain wind speed from from IEC files
else
   VW = FE.wind{2,1}(1:3)';  % Wind Speed Vector from *WIND card
end
VT = VW - VB ;   % Resultant Velocity seen by the blade
Va =  VT' * en;  % Total axial velocity
Vt =  VT' * et;  % Total tangential velocity
VWn = VW' * en;  % Normal to rotor wind velocity



aax=0.0;     atg=0.0;    ptlf=1;   nb=3; % Induction factors, Prandtl Tip Loss Factor and Fixed Number of Blades  
%------------------------------------------------------------------------
%                      ACTIVATE BEM ITERATION
%------------------------------------------------------------------------
if strcmp(FE.rotor{1,3},'YES')
    niters=20;
else
    niters=1;
end


%------------------------------------------------------------------------
%                         ROTOR STARTING
%------------------------------------------------------------------------
rotorstart =0; 
if FE.rotor{2,2}(1) == 1  % if start flag is active
    if X.tstep <= FE.rotor{2,2}(2)  % if time step is between start steps
        rotorstart = 1;
        niters = 1;
        startforce = FE.rotor{2,2}(3);
    end 
end



%------------------------------------------------------------------------
%                         BEM ALGORITHM
%------------------------------------------------------------------------

chord = interp1( btb(:,1) , btb(:,4) , rm, 'linear' );

for i=1:niters;

    % VER PORQUE PUSE kax y ktg ADENTRO DEL LOOP DE ITERACIÓN SI NO
    % ESNECESARIO...
    kax = 2*pi*(re^2-ri^2)*rho*VWn^2 * ptlf /nb   ; 
% 	T1D = kax * aax * (1-aax);  % 1D Momentum Thrust per blade element
    
    ktg = pi*(re^4-ri^4)*rho*ome*VWn^2  * ptlf /nb;
% 	M1D = ktg * atg * (1-aax);  % 1D Momentum Torque  per blade element
    
    Vi = Va*en *(1-aax) + Vt*et * (1+atg);  % Global Incidence Velocity Vector 
    
    vi = ei' * Vi;    % Local Incidence Velocity Vector 
    
    mvi = sqrt ( vi(2)^2 + vi(3)^2 ); % Squared Modulus of the in-plane velocity
    
    alpha = atan2(vi(3),vi(2)) * 180/3.14159;  % Four Cuadrant Angle of attack in degrees 
    
    acoef = interp1( ftb(:,1) , ftb(:,2:3) , alpha, 'linear' );   % Spline Interpolated CL 
    CL = acoef(1);    CD = acoef(2);
    
    uforce = 0.5 * rho * chord * mvi^2;
    LL =  CL * uforce; % Local Lift 
    DL =  CD * uforce; % Local Drag
    
    FC = [0 ; DL ; LL ];         % Force Vector for Lift and Drag (rotated alpha)
    Ra = expmap( [alpha*pi/180 0 0] ); % Rotation matrix for alpha
    FAL = Ra * FC; % Local Aerodinamic force per unit length
     
    % STARTING FORCE, CHECK should be outside the routine
    if rotorstart == 1
        FAL = [0 startforce 0]';
    end
    
    FAG = ei * FAL;    % Global Aerodinamic force per unit length
    T2D =  dr * FAG' * en; % Bidimensional thrust in the element
    M2D = f_cross( xcr, dr*FAG )' * en;
    
%     aaxnew = ( -1 + sqrt(1-4* ( abs(T2D)/k1 ) ) ) / -2;        
    % FIXED POINT ITERATION for induction factors
    aaxnew = aax^2 + abs(T2D) / kax;  %  Slower form: aax = abs(T2D) / ( k1 * (1 - aax) )
    atgnew = aax*atg+ abs(M2D / ktg);
    
 
    % CONVERGENCE TEST    
    if i>=2
        atol=1E-3;
        if  abs(aaxnew - aax)  <= atol && abs(atgnew - atg)  <= atol % Convergence criterion
            aax = aaxnew;
            atg = atgnew;
        break
        end
    end 
    
    aax = aaxnew;
    atg = atgnew;
    
        if i>=2 && i>=niters
             disp('FORCING INDUCTION FACTORS CONVERGENCE')
        end
end


fa = [ FAG ; 0 ; 0 ; 0 ];     % Element Aerodynamic force vector per unit length



X.aero(el,1:7,X.tstep+1) = [ome aax atg i T2D M2D alpha];  % Info table

%---------------------------------------------------------------
%                        WARNING FLAGS
%---------------------------------------------------------------
if X.wflag(1)==0         
     disp ('AEROWARNING: Rotor Normal Vector is Fixed')
     disp (['AEROWARNING: BEM Geometry is set to ' metodo])
     disp ('AEROWARNING: Yaw angle is set to zero in iecspeed.m ')
     disp ('AEROWARNING: Update model for hub height different from zero in iecspeed.m ')
     X.wflag(1)=1;
end



