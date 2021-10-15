%============================================================
%  Purpose: arc length method
%
%  Martin Saravia - 20-03-2011
%  Note: arclength 2 replaces the old arclength 1, vectors and time stepping info is
%        slightly different than arclength 1 introduced in nlbeam. This modified function 
%        was introduced in mbdyn1.2
%============================================================ 

function [X,DU]=arclength2(X,K,DU,U,F)


% Desplazamiento por CARGA TOTAL con KT en convergencia del paso anterior
% UQ=K\Fr;        
[UQ(X.adof,1),X]=ldlsolve(X,K,F(X.adof,1),3);

% Desplazamiento total correspondiente al residual
if X.iter==1 
      DUT=DU;
    else

    % Rutina vieja donde cambio la dimension de los vectores por bc de eliminacion, ahora no lo uso mas
    % directemento asigno los dof que quiero del vector total con el vector de grados de libertad activos
    % X.adof.
    %   if X.tstep==1   
    %     DUT=bcelifwd(U(:,X.tstep+1),X);
    %   else
    %     DUT=U(:,X.tstep+1)-U(:,X.tstep); %OJOOOOOO -- NO DEBO SUMAR DU, PORQUE DUT DE LA ITERACION ANTERIOR
    %   end
    %   DUT=bcelifwd(DUT,X); 

    DUT=U(:,X.tstep+1)-U(:,X.tstep);
end


%%-------------------------------------------------------------------------------
%                                  PREDICTOR 
%--------------------------------------------------------------------------------
if X.iter==1
  X.dlambda = X.predsign * X.arclength / (sqrt(UQ'*UQ));  %ojo, incluyo las rotaciones
  DU = DU + X.dlambda*UQ;                                 %ojo, no va UQ disp, sino el vector real.
  X.lambda = X.T(3,X.tstep-1) + X.dlambda;
  X.T(3,X.tstep) = X.lambda; 
end

%%-------------------------------------------------------------------------------
%                                 CORRECTOR
%--------------------------------------------------------------------------------
if X.iter>=2
  aa = UQ'*UQ;
  bb = 2*UQ'*(DUT+DU);
  cc = ( (DUT + DU)'*(DUT + DU) ) - X.arclength^2;
  dlambda1 = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
  dlambda2 = (-bb-sqrt(bb^2-4*aa*cc))/(2*aa);
  X.dlambda1=dlambda1 ;
  X.dlambda2=dlambda2 ;
  
  if bb^2-4*aa*cc<=0
    X.stop=1;
    disp('NO HAY RAICES REALES PARA EL PREDICTOR')
  end
  
    % Crisfield corrector root selection    
    %  a4 = DUTd'*DUd + DUTd'*DUTd;   
    %  a5 = DUTd' * UQd;
    %   if (a4+a5*dlambda2) >= (a4+a5*dlambda1) % Crisfield
    %     X.dlambdac=dlambda2;
    %     else
    %     X.dlambdac=dlambda1;
    %   end 
      
    % Ritto-Correa corrector root selection
    tt=DUT' * UQ ;                    
    if (tt*dlambda1) <= (tt*dlambda2)
    X.dlambda=dlambda2;
    else
    X.dlambda=dlambda1;
    end


  DU = DU + X.dlambda * UQ;
  X.lambda = X.lambda + X.dlambda;
  X.T(3,X.tstep) = X.lambda; 
end


%--------------------------------------------------------------------------------
%    DETERMINACIÓN DE LA DIRECCIÓN PARA EL PRÓXIMO PREDICTOR
%--------------------------------------------------------------------------------
if X.iter==1
  X.DAQ=0; %inicializo 
end

% Dirección sin rotaciones
X.DAQ=(X.dlambda*UQ)+X.DAQ;    %calculo el desp total correspondiente a dlambda
DAR=DUT;
direction=UQ'*(X.DAQ+DAR);     %Predictor de signo del predictor

% %Dirección con rotaciones
% X.DAQ=(X.dlambda*UQ)+X.DAQ;   %calculo el desp total correspondiente a dlambda
% DAR=DUT;
% direction=UQ'*(X.DAQ+DAR);    %Predictor de signo del predictor

if direction>=0 
  X.predsign=1;
else
  X.predsign=-1;
end 

X.T(4,X.tstep)=X.predsign;