%%-------------------------------------------------------------
%
% U: Total Displacements
% DU=Di: Incremental Displacements
% dU: Iterative displacements
%
%---------------------------------------------------------------

function [X,U,Ui,Vel,Acc,Ome,Gam,RT] = iterate (FE,X,N,U,Ui,Vel,Acc,Ome,Gam,K,C,M,F,FD,FI,FM,FA,RT)

dU=zeros(X.sdof,1);
%============================================================
%                      ESTATICA
%============================================================
if strcmp(FE.step{1,3},'STATIC')==1  

    %  -------------    RESIDUAL --------------
    disp (norm(FA(X.adof,X.tstep+1)))
    disp (X.tstep)
    Res =  F(X.adof,X.tstep+1) + FD(X.adof,X.tstep+1) + FA(X.adof,X.tstep+1) - FI(X.adof,X.tstep+1); 

    dU(X.adof,1) = K \ Res;
    
    %   ------------  ARCLENGTH --------------
    if X.ArcL==1
       [X,dU]=arclength2(X,K,dU,U,F);
    end

    %  ----------------  U update -------------------
    Ui(X.adof,X.tstep+1) = Ui(X.adof,X.tstep+1) + dU(X.adof,1);       % Actualizacion iterativa 
    U(X.adof,X.tstep+1) = U(X.adof,X.tstep) + Ui(X.adof,X.tstep+1);

    %   ----------  CONVERGENCE CHECK --------------------
    [X] = convcheck(X,FE,Res,F,FD,FA,dU,U);  % Residual de la configuracion (sin LMs)

    RT = rotupds(X,RT,Ui,'corr');
end    


%============================================================
%                      DINAMICA
%============================================================   
if strcmp(FE.step{1,3},'DYNAMIC')==1  
    
    method = FE.step{2,2}(1);   dt = FE.time(3);   rho = FE.step{2,2}(2);
    
    %===============================================================================
    %                               RESIDUAL    
    %===============================================================================
     if  FE.step{2,2}(3) == 1   % NLGEOM=YES
        Res = F(X.adof,X.tstep+1) + FD(X.adof,X.tstep+1) + FA(X.adof,X.tstep+1) - FM(X.adof,X.tstep+1) - C*Vel(X.adof,X.tstep+1) - FI(X.adof,X.tstep+1);
     end
     if FE.step{2,2}(3) == 0   % NLGEOM=NO
        Res = F(X.adof,X.tstep+1) + FD(X.adof,X.tstep+1) + FA(X.adof,X.tstep+1) - M*Acc(X.adof,X.tstep+1) - C*Vel(X.adof,X.tstep+1) - K*U(X.adof,X.tstep+1);
     end
     

     
    %===============================================================================
    %                                  ITERATION 
    %===============================================================================    


    %------------- INTEGRATION METHOD ----------------------------
    if method==1 || method==2 % ALPHA 
        alpham=(2*rho-1)/(rho+1);      alphaf=rho/(rho+1);
        gamma=0.5+alphaf-alpham;       beta=0.25*(gamma+0.5)^2; 

        gammap=gamma/(dt*beta);        
        betap=(1-alpham)/(dt^2*beta*(1-alphaf));
    end
    
    if method==4  % NEWMARK UL 
        gammap=(2/dt);        
        betap=(4/dt^2);
    end   

    %------------- ITERATION ----------------------------
    S =  betap*M + gammap*C + K ;     % ITERATION MATRIX   
    dU(X.adof,1) = S \ Res;           % SOLUTION OF INCREMENTAL DISPLACEMENT

    % ------------- DYNAMIC VARIABLES UPDATE ------------
    Ui(:,X.tstep+1) = Ui(:,X.tstep+1) + dU;  % Aca actualizo los multiplicadores
    U(:,X.tstep+1) = U(:,X.tstep) + Ui(:,X.tstep+1);
    
    Vel(X.sfdof,X.tstep+1) = Vel(X.sfdof,X.tstep+1) + gammap * dU(X.sfdof);
    % Vel(FE.ic(1,:),X.tstep+1) = FE.ic(2,:) ; %Velocidad impuesta;
    Acc(X.sfdof,X.tstep+1) = Acc(X.sfdof,X.tstep+1) + betap * dU(X.sfdof);
    
    [X] = convcheck(X,FE,Res,F,FD,FA,dU,U);  % Residual de la configuracion (sin LMs)
    
    [RT,Ome,Gam] = rotupdd(X,RT,Ui,Vel,Acc,Ome,Gam,'corr');   % Actualización de rotaciones y velocidades angulares
           

    
end 


