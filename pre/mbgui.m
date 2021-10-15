function PH = mbgui (FE)

%==========================================================================
%                          INICIALIZACION GRAFICA
%==========================================================================
 
% save([FE.Nombre '.fem'],'FE','N','E','C','BC','Fs')       % Guardo el trabajo.


PH.out=zeros(2,1); %inicializo la tabla de convergencia
 
PH.f1=figure('Position',[10 50 1300 600], 'Name','GRAFICAS');
%textos de informacion para rutina messages


% Tabla de convergencia
% cnames = {'step','LPF','iters'};
% uitable('Parent',f2,'Data',[],'ColumnName',cnames,'Units','Normalized','Position',[0.05 0.05 0.95 0.45],'Tag','tabla');


% Figura de informacion de convergencia
% f2=figure('Position',[930 100 400 500], 'Name','INFORMACION');

% Plot Figures
PH.uip2 = uipanel('Parent',PH.f1,'Title','Real Time Plots','FontSize',10,'BackgroundColor','white','Position',[0.01 0.01 0.78 0.98]);

% Run Information Panel
PH.nfopnl = uipanel('Parent',PH.f1,'Title','Run Information Panel','FontSize',10,'BackgroundColor','white','Position',[0.8 0.01 0.19 0.68]);



PH.uic1=uicontrol('Parent',PH.nfopnl,'Style', 'text', 'String', ' ','Position', [10 40 210 20]); % STEP - ITERACION
PH.uic2=uicontrol('Parent',PH.nfopnl,'Style', 'text', 'String', ' ','Position', [10 10 210 20]); %CONVERGENCIA
PH.uic3=uicontrol('Parent',PH.nfopnl,'Style', 'text', 'String', ' ','Position', [10 70 210 20]); %DATOS
PH.axit=axes('Parent',PH.nfopnl,'Position', [ 0.1 0.2 0.88 0.78]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











