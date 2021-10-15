%-------------------------------------------------
%      GRAFICADOR DE CURVAS Y EXPORTADOR A EPS
%
% Martin Saravia - 28-04-2010
%
% Nota: Requiere GhostScript, export_file, pdftops
%       Ver export_file
%-------------------------------------------------
clc
clear all
addpath([cd '/pre/']); config ; 
addpath([cd '/plot/']);
%-------------------------------------------------
%                  PARAMETROS
%-------------------------------------------------
F1=1;
%-------------------------------------------------
%                MATLAB DATASET NLBEAM 1.6 (euleriana)
%-------------------------------------------------
load('bipendulo_2m_20el_EGlass.mat')
tstep=400;
umatE=[ 0 U(129,1:tstep/2)]; 
% vmatE=[ 0 -U(302,:)];
% wmatE=[ 0 -U(303,:)];
lfmatE=[ 0 (X.T(2,1:tstep/2))];

%-------------------------------------------------
%                MATLAB DATASET NLBEAM 2.0 (Total lagrangiana)
%-------------------------------------------------
% load('arc45_lam45_TL.mat')
% 
% umatTL=[ 0 -U(301,:)]; 
% vmatTL=[ 0 -U(302,:)]; 
% wmatTL=[ 0 -U(303,:)]; 
% lfmatTL=[ 0 F1*(X.T(3,:))];



%-------------------------------------------------
%              ABAQUS DATASET
%-------------------------------------------------

ABA=load('bipenduloEGlass.rpt','');
uaba=ABA(1:tstep,2); % promedio de los desp en los ptos elegidos
% vaba=-( (ABA(:,6)+ABA(:,7)+ABA(:,8)+ABA(:,9))/4 );
% waba=-( (ABA(:,10)+ABA(:,11)+ABA(:,12)+ABA(:,13))/4 );
lfaba=ABA(1:tstep,1);


h=figure(1);
set(h,'Position', [200, 200, 500, 300])
set(h, 'color', 'white')

hold on

p1=plot(lfaba,uaba);
set(p1, 'LineStyle', '-.', 'LineWidth', 0.5, 'Color', 'Black');
% set(p2, 'Marker', 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 1.0);

p2=plot(lfmatE,umatE);
set(p2, 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'Black');
% set(p1, 'Marker', 'x', 'MarkerFaceColor', [0.5 0.5 0.5],
% 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 1.4);



set(gca,'FontSize',8,'FontName','Times','Box','on')
xlabel('Time (s)','FontSize',8,'FontName','Times'); 
ylabel('Displacement','FontSize',8,'FontName','Times'); 
leg = legend('Abaqus Shell','Present',2);
set(leg,'Interpreter','none','Location','NorthWest')

%-------------------------------------------------
%             EXPORTACION FIGURA
%-------------------------------------------------
export_fig ('bipenduloEGlassh.eps',-eps, h) 
