function [X]=messages(handles,X,FE)


% % Tabla de Convergencia
% X.CT(X.tstep,1:3)=[X.tstep X.T(3,X.tstep) X.iter];
% data_tabla=num2cell(X.CT);
% find_tabla=findobj('Tag', 'tabla');
% set(find_tabla, 'Data', data_tabla); 

% Convergence Information

set(handles.sirtext,'String',[num2str(X.tstep) '                ' num2str(X.iter)   '               '    '0']);

% msg1=['Step: ' num2str(X.tstep)];
% msg2=['Iteración: ' num2str(X.iter)];
% set(PH.uic1, 'String', [msg1 '            ' msg2]);
% 
% 
% % Rotor Information
% if isfield (FE,'rotor') ==1  
%     el = FE.rotor{2,1}(2);
%     msg3=['Omega: ' num2str(X.aero(el,8,X.tstep+1))];
%     msg4=['Power: ' num2str(X.aero(el,9,X.tstep+1))];
%     find_uit3=findobj('Tag', 'uit2');
%     set(find_uit3, 'String', [msg3 '            ' msg4]);   
% end
% 
% 
% if X.conv==1
%     set(PH.uic2, 'String', ' CONVERGENCE');
% else
%     set(PH.uic2, 'String', ' NO CONVERGENCE!');
% end
% 
% 
plot(handles.axes5,X.stpnfo(:,1), X.stpnfo(:,3) )



