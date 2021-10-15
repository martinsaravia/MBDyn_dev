%% ========================================================================
%                      TRIADS FOR WIND TURBINE BLADE
% ==========================================================================

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  BE CAREFULL ONLY WORKS BEFORE SECTION ASSIGNMENTS, didx= E(el,3) is only
%  valid before that assigments.

function [A] = triadsb(FE,N,E,el)

%--------------------------------------------------------------------------
nidx = E(el,8)*2-1;  
rootnode = FE.rotor{2,1}(1);  bladeshape = FE.sec{nidx,6};
%--------------------------------------------------------------------------
nn1=E(el,4);  nn2=E(el,5);  nn3=E(el,6);  % Numero de nodo del elemento
xn1 = N(nn1,2:4);  xn2 = N(nn2,2:4);  xn3 = N(nn3,2:4);    % Node Positions

xroot = N(rootnode,2:4); 
%--------------------------------------------------------------------------

xc = 0.5*xn1 + 0.5*xn2;                      % Element Center Position
rpos = norm(xc-xroot);    % Position relative to root

% Rotated third node coordinates according to pitch
bldaxs = (xn2 - xn1)/norm(xn2-xn1); % blade unit axis

pitch = interp1( bladeshape(:,1) , bladeshape(:,2) , rpos,'cubic' );   % Spline Interpolated Pitch
pitchR = expmap ( (pitch*pi/180)*bldaxs ); % Pitch exponential map NOTE THE plus SIGN IN ORDER TO ROTATE THE FOIL ANTICLOCKWISE AND GIVE POSITIVE ATTACK
xn3 = pitchR * xn3'; % new coordinates



% DIRECTORS

d1 = (xn2-xn1)'/norm(xn2-xn1); % Director in direction 1

dt = (xn3'-xn1)';  %Temporary director (to define the plane)

d2 = f_cross(dt,d1);  d2=d2/norm(d2);

d3 = f_cross(d1,d2);  d3=d3/norm(d3);   % 

A = [d1 d2 d3];




