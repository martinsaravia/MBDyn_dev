%-------------------------------------------------
%      GRAFICADOR DE CURVAS PARA PAPER
%
% Martin Saravia - 07/2012
%
% Nota: Requiere GhostScript, export_file, pdftops
%       Ver export_file
%-------------------------------------------------
function plot2paper(var1,var2,filename,export,color,style,namex,namey)

h=figure(1);
set(h,'Position', [200, 200, 600, 300])
set(h, 'color', 'white')


p1=plot(var1,var2);
set(p1, 'LineStyle', style, 'LineWidth', 1.0, 'Color', color);
% set(p1, 'Marker', 'x', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 1.4);

set(gca,'FontSize',10,'FontName','Times','Box','on')
xlabel(namex,'FontSize',10,'FontName','Times'); 
ylabel(namey,'FontSize',10,'FontName','Times'); 

if strcmp(export,'YES')
 export_fig ( (strcat(filename,'.eps')),-eps, h) 
end
