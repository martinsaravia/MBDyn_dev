%============================================================
%  PURPOSE: Solver linear system
%  IN: Matrix and vector and type
%  OUT: Vector and diagnosis
%
%  Martin Saravia - 12-04-2010(1.5) -
%
% 
%============================================================

function [DU,X]=ldlsolve(X,K,vec,type)


%===========================================================
%              TIPO 1:  DESCOMPOSICION LDL
%===========================================================
if type==1
    [Lm, Dm, pm] = ldl(K, 'vector');

  %   fprintf(1, 'The error norm ||M(pm,pm) - Lm*Dm*Lm''|| is %g\n', ...
  %     norm(K(pm,pm) - Lm*Dm*Lm'));

    DU(pm,:) = Lm'\(Dm\(Lm\(vec(pm,:))));

  %   fprintf('The absolute error norm ||x - ones(size(b))|| is %g\n', ...
  %     norm(DU - ones(size(vec))));

  % NEGATIVE PIVOTS DETERMINATION
%   for i=1:X.rdof
%     pivot=Dm(i,i);
%     pivot2=min(min(Dm));
%     if pivot<=-1e-15
%       X.negpivot=-1;
%       fprintf('I found a negative pivot: Magnitude is %g\n',pivot)
%     end
%   %   if pivot2<=-1e-15
%   %     X.negpivot=-1;
%   %     fprintf('I found a negative off diagonal pivot: Magnitude is %g\n',pivot2)
%   %   end
%   end
% 

 end

%===========================================================
%              TIPO 2:  DESCOMPOSICION LU
%===========================================================
if type==2  
  [Lf,Uf] = lu(sparse(K));  % Factorizo la K tangente
  DU = Uf\(Lf\vec);   % Desp incremental
end


%===========================================================
%              TIPO 3:  RESOLUCION DIRECTA CON OPERADOR \
%===========================================================
if type==3 
  DU=K\vec;        %% Desp Incremental calculado con back slash
end


%===========================================================
%              TIPO 4:  LINSOLVE
%===========================================================
if type==4
  DU=linsolve(K,vec);   %% Resuelvo con LU pero con rutina de Matlab
end



% options.disp=0; vals=5;
% ev=eigs(K,vals,'SM',options);
% X.eigs(1:vals,X.tstep)=ev;
% X.eigs(vals+1,X.tstep)=X.lambda;

