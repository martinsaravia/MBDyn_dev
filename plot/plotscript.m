
% 

% %------------------------------------------------------------
% %                    DEFORMED SHAPES
% % %------------------------------------------------------------

% clc
% clear all
% fstep=1000;
% load('SNLWT100A2_NWP113_Lin')
% scale=[-100 100];
% % spatch1d(X,E,N,U,1,{'PATCH' 353},{'' ''},gscale,view)
% % spatch1d(X,E,N,U,1,{'PATCH' 1},{'' ''},gscale,view)
% view=([0 90]); 
% drawnow



% %------------------------------------------------------------
% %                    VARIAS
% % %------------------------------------------------------------
% clc
% clear all
% fstep=1000;
% load('SNLWT100A2_NWP113')
% var1=X.T(2,1:fstep); 
% % var2=squeeze(SE(2,1,1:fstep)); % Bending Moments
% % var2=squeeze(X.aero(64,9,1:fstep))/1E6; % Bending Moments
% var2=U(308,1:fstep);% Tips Displacement
% % var2=squeeze(X.aero(33,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black','-','','');    %axis([0 30 0 25])
% hold on 

% export_fig ( (strcat('SNLWT100A2_NWP113_nodamp_defshape','.eps')),-eps, figure(1))


%------------------------------------------------------------
%              COMPARACION DAMP - NODAMP
%------------------------------------------------------------
% clear all
% fstep=1000;
% load('SNLWT100A2_EWM50_damp')
% var1=X.T(2,1:fstep); 
% % var2=squeeze(SE(2,1,1:fstep)); % Bending Moments
% % var2=squeeze(X.aero(64,9,1:fstep))/1E6; % Bending Moments
% var2=U(312,1:fstep)*180/pi;% Tips Displacement
% % var2=squeeze(X.aero(33,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black','-','','');    %axis([0 30 0 25])
% hold on  


% 
% clear all
% fstep=1000; namex='Tiempo (s)'; namey='Giro Torsional (grados)';
% load('SNLWT100A2_EWM50_damp_Lin')
% var1=X.T(2,1:fstep); 
% % var2=squeeze(SE(2,1,1:fstep));% Bending Moments
% % var2=squeeze(X.aero(64,9,1:fstep))/1E6; % Bending Moments
% var2=U(312,1:fstep)*180/pi;% Tips Displacement
% % var2=squeeze(X.aero(33,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black',':',namex,namey);    %axis([0 30 0 25])
% 
% export_fig ( (strcat('SNLWT100A2_EWM50_ThetaTorsTip_linvsnolin','.eps')),-eps, figure(1))


%------------------------------------------------------------
%              COMPARACION ENTRE PALAS
%------------------------------------------------------------
% %  ANGLE OF ATTACK PLOT
clc
clear all



load('SNL100/results/SNLWT100A2_NWP113_nodamp')
fstep=2:500;
var1=X.T(2,fstep); 
var2=squeeze(SE(3,[1 3 5],fstep)); % Bending Moments
% var2=U(302,1:fstep);% Tips Displacement
% var2=FA(129,1:fstep);% Tips Displacement
% var2=FI(6,1:fstep);% Tips Displacement
% var2=squeeze(X.aero(38,7,1:fstep)); % Angle of attack
% var2 = squeeze ( EE(2,1,fstep) );
plot2paper(var1,var2,'NN','NO','Red','-','','');    %axis([0 30 0 25])
hold on 

clear all

load('SNL100/results/SNLWT100A2_NWP113_Lin_nodof4_nodamp_S6new')
fstep=2:500;
var1=X.T(2,fstep); 
var2=squeeze(SE(3,[6 8 10],fstep)); % Bending Moments
% var2=U(302,1:fstep);% Tips Displacement
% var2=FA(129,1:fstep);% Tips Displacement
% var2=FI(6,1:fstep);% Tips Displacement
% var2 = squeeze ( EE(2,1,fstep) );
% var2=squeeze(X.aero(38,7,1:fstep)); % Angle of attack
plot2paper(var1,var2,'NN','NO','Blue','-','','');    %axis([0 30 0 25])
hold on 

% clear all
% fstep=300;
% load('run/linearSCBLT')
% 
% var1=X.T(2,1:fstep); 
% % var2=squeeze(SE(6,54,1:fstep)); % Bending Moments
% var2=U(69,1:fstep);% Tips Displacement
% % var2=squeeze(X.aero(17,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black','-','','');    %axis([0 30 0 25])
% hold on 


% % clear all
% % load('SNLWT100A2_NWP113_nodamp')
% var1=X.T(2,1:fstep); 
% % var2=squeeze(SE(2,3,1:fstep));% Bending Moments
% % var2=Acc(308,1:fstep);% Tips Displacement
% var2=squeeze(X.aero(33,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black','--');    %axis([0 30 0 25])
% 
% var1=X.T(2,1:fstep); 
% var2=squeeze(SE(2,5,1:fstep)); % Bending Moments
% % var2=Acc(302,1:fstep); % Tips Displacement
% var2=squeeze(X.aero(51,7,1:fstep)); % Angle of attack
% plot2paper(var1,var2,'NN','NO','Black',':');    %axis([0 30 0 25])


% export_fig ( (strcat('SNLWT100A2_Extreme25ms','.eps')),-eps, figure(1))

 



