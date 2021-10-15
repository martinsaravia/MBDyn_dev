%============================================================
%  Purpose: convergence check
%  IN: residual and external force
%  OUT: convergence flag and iteration log
%
%  Martin Saravia - 29-03-2010
%============================================================

function [X]=convcheck(X,FE,Res,F,FD,FA,dU,U)

if FE.step{2,2}(3) == 0   % NLGEOM=NO
 X.conv = 1;
 return
end

ftol=1E-4; % Tolerancia 
utol=ftol;

resnorm=sqrt(Res'*Res);    % Euclidean Norm del residual

if strcmp(FE.step{1,3},'STATIC') == 1
   fext = F(X.adof,X.tstep+1) + FD(X.adof,X.tstep+1) ;  % Static Residual
else
   fext = F(X.adof,X.tstep+1) + FD(X.adof,X.tstep+1) + FA(X.adof,X.tstep+1);   % Dynamic Residual
end
fextnorm=sqrt( fext' * fext);  % Euclidean Norm de la fuerza externa

dUnorm=sqrt(dU'*dU);
Unorm=sqrt(U(X.adof,X.tstep+1)'*U(X.adof,X.tstep+1));


if X.iter>=FE.param(5)  %obligo a iterar ciertas veces
    if  dUnorm<=(utol*Unorm) ||  resnorm<=(ftol*fextnorm)   
        X.conv=1;   
    end
end

if X.iter>=FE.param(4)
    X.conv=1; 
    fprintf(fopen(X.outputfile,'a'),'WARNING ! Forcing convergence in time step %d. \n\n', X.tstep);
    fclose(fopen(X.outputfile,'a'));
end


if X.iter>=FE.param(4)
    if X.conv==0
    fprintf(fopen(X.outputfile,'a'),'ERROR !!!  No convergence at step: %d \n\n', X.tstep);
    error('xxxxx  ERROR !!!  No convergence !!!  xxxxx')
    end
end

X.stpnfo(X.tstep+1,:) = [X.tstep  X.iter  resnorm   dUnorm ];

