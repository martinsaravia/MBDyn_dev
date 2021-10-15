function [X,DU]=arclength(X,K,DU,U,F)


% Desplazamiento por carga total con la tangente en convergencia del paso anterior
% UQ=K\Fr;        
[UQ,X]=ldlsolve(X,K,F(X.adof,1),3);

% Desplazamiento total correspondiente al residual
if X.iter==1 
  DUT=DU;
else
    
%   if X.tstep==1   
%     DUT=bcelifwd(U(:,X.tstep+1),X);
%   else
%     DUT=U(:,X.tstep+1)-U(:,X.tstep); %OJOOOOOO -- NO DEBO SUMAR DU, PORQUE DUT DE LA ITERACION ANTERIOR
%   end
%   DUT=bcelifwd(DUT,X); 
  
DUT=U(:,X.tstep+1)-U(:,X.tstep);
  
  
  
  
end

% Elimino las rotaciones de todos los vectores
% UQd=bceliback(UQ,X);   
% DUTd=bceliback(DUT,X);   
% DUd=bceliback(DU,X);
% for i=X.nds        
%   rdof1=i*6-2; rdof3=i*6;
%   UQd(rdof1:rdof3,1)=0;
%   DUTd(rdof1:rdof3,1)=0;
%   DUd(rdof1:rdof3,1)=0;
% end
% UQd=bcelifwd(UQd,X);   
% DUTd=bcelifwd(DUTd,X);   
% DUd=bcelifwd(DUd,X);

%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                                  PREDICTOR 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if X.iter==1
  X.dlambda = X.predsign * X.arclength / (sqrt(UQd'*UQd));  %ojo, no incluyo las rotaciones
  DU = DU + X.dlambda*UQ; %ojo, no va UQ disp, sino el vector real.
  X.lambda = X.T(3,X.tstep-1) + X.dlambda;
  X.T(3,X.tstep) = X.lambda; 
end

%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                                 CORRECTOR
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if X.iter>=2
  aa = UQd'*UQd;
  bb = 2*UQd'*(DUTd+DUd);
  cc = ( (DUTd + DUd)'*(DUTd + DUd) ) - X.arclength^2;
  dlambda1 = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
  dlambda2 = (-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
  X.dlambda1=dlambda1 ;
  X.dlambda2=dlambda2 ;
  
  if bb^2-4*aa*cc<=0
    X.stop=1;
    disp('NO HAY RAICES REALES PARA EL PREDICTOR')
  end
    
%  a4 = DUTd'*DUd + DUTd'*DUTd;    % Crisfield corrector root selection
%  a5 = DUTd' * UQd;
%   if (a4+a5*dlambda2) >= (a4+a5*dlambda1) % Crisfield
%     X.dlambdac=dlambda2;
%     else
%     X.dlambdac=dlambda1;
%   end 
      
  
  tt=DUTd' * UQd ;                     % Ritto-Correa corrector root selection
  if (tt*dlambda1) <= (tt*dlambda2)
    X.dlambda=dlambda2;
    else
    X.dlambda=dlambda1;
  end




  DU = DU + X.dlambda * UQ;
  X.lambda = X.lambda + X.dlambda;
  X.T(3,X.tstep) = X.lambda; 
end



%   X.DU146=DUd(146,1);
%   X.UQ146=UQd(146,1);
%   X.DUT146=DUTd(146,1);
%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                                 ACTUALIZACIÓN
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                  


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% DETERMINACIÓN DE LA DIRECCIÓN PARA EL PRÓXIMO PREDICTOR
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if X.iter==1
  X.DAQ=0; %inicializo 
end

% Dirección sin rotaciones

X.DAQ=(X.dlambda*UQd)+X.DAQ;    %calculo el desp total correspondiente a dlambda
DAR=DUTd;
direction=UQd'*(X.DAQ+DAR); %Predictor de signo del predictor

% %Dirección con rotaciones
% X.DAQ=(X.dlambda*UQ)+X.DAQ;    %calculo el desp total correspondiente a dlambda
% DAR=DUT;
% direction=UQ'*(X.DAQ+DAR); %Predictor de signo del predictor



if direction>=0 
  X.predsign=1;
else
  X.predsign=-1;
end 

X.T(4,X.tstep)=X.predsign;