%-----------------------------------------------------
%     CALCULATE WIND SPEED FROM IEC FILES
%
% WARNING ! Routine only valid for 0 yaw angle.
% Hubnode heigh is variable, good to note
%-----------------------------------------------------

function VW = iecspeed (FE,X,wdata,xcr,xhub,Ea,Ev,Et)


yaw =0;

RR = FE.rotor{2,1}(9); % rotor Radius
hheight = FE.rotor{2,1}(10); % rotor Radius

Vhub = wdata(2);
wdir = wdata(3);
Vshr = wdata(4);
Hshr = wdata(5);
Vshrexp = wdata(6);
VshrL = wdata(7);
Vgust = wdata(8);

vang = wdir-yaw;

xvb = Ev * xcr + hheight;  % Vertical Coordinate of the Blade
xtb = Et * xcr;  % Horizontal Coordinate of the Blade
xab = Ea * xcr;  % Relative Axial position of the Blade
xvh = Ev * xhub + hheight; % Vertical Coordinate of the Hub



V1 = Vhub * ( (1+xvb/xvh)^Vshrexp + Hshr*( xtb*cosd(vang)+(xab*sind(vang))/(2*RR) + VshrL*xvb)) + Vgust;

Vx = V1 * cosd (vang);
Vy = V1 * sind (vang);

VW = (-Vx*Ea + Vy*Et)';



