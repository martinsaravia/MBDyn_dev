%============================================================
%  PURPOSE: Update Everything
%  IN: vector
%  OUT: updated vectors
%
%  Martin Saravia - 10/08/2010
%
%============================================================

function    [F,X,U,Ui,Vel,Acc,Ome,Gam,A0,RT]=advance(X,FE,F,U,Ui,Vel,Acc,Ome,Gam,A0,RT)

X.conv=0; 

% ------------------------------------------------------------------------
%                            S T A T I C S
% ------------------------------------------------------------------------

if strcmp(FE.step{1,3},'STATIC')==1

    %   TIME EVOLUTION INFORMATION
    X.T(1,X.tstep+1) = X.tstep; 
    
    % ARCLENGTH
    if strcmp(FE.step{2},'arclth')==1
      if X.tstep>=FE.dyn{4}
        X.flag(11)=1;
        X.T(3,X.tstep) = X.lambda;                   % lambda
      end
    end
    
    if X.flag(11)==0
        X.dlambda=FE.time(3);
        X.lambda=X.dlambda*X.tstep;
    end

    X.T(2,X.tstep+1) = X.lambda*FE.time(2);            % t actual
    X.T(3,X.tstep+1) = X.lambda;                       % lambda

    % INREMENTAL UPDATE OF CONFIGURATION VARIABLES (INCLUDING IMPOSED CONFIGURATIONS)
    Ui(X.elidof,X.tstep+1) = FE.bc (2,:)';                                       % Impose incremental displacements
    U(X.elidof,X.tstep+1) = U(X.elidof,X.tstep) + Ui(X.elidof,X.tstep+1);        % Impose incremental displacements  
    U(X.adof,X.tstep+1) = U(X.adof,X.tstep);                                     % Actualizacion de desp
    Ui(X.adof,X.tstep+1) = 0;
    
    RT = rotupds(X,RT,Ui,'pred');  % Rotation update for the prediction
end


% ------------------------------------------------------------------------
%                         P R E D I C T O R -  D Y N A M I C S
% ------------------------------------------------------------------------
if strcmp(FE.step{1,3},'DYNAMIC')==1
    
    method = FE.step{2,2}(1); 
    dt = FE.time(3);  

%     % VELOCIDADES ANGULARES INICIALES
%     for ni = 1: length ( X.actnodes )
%         nn = X.actnodes (ni);
%         idxn = nn*3-2 : nn*3;
%         rdofn = X.CTable(nn,4:6);
%         Tn  =  RT(idxn,1:3,2);                        % Incremental Tangential Map
%         Acc(rdofn,X.tstep) = Tn * Acc(rdofn,X.tstep) + xidT(Vel(rdofn,X.tstep),Ui(rdofn,X.tstep)) * Vel(rdofn,X.tstep);
%         Vel(rdofn,X.tstep) = Tn * Vel(rdofn,X.tstep); % Angular Velocity Update
%     end
    
    
    %------------------------------------------------------------------------------
    %                                 ALPHA DE CARDONA
    %------------------------------------------------------------------------------
    if method == 1; %ALPHA METHOD CARDONA 

        rho = FE.step{2,2}(2);
        alpham=(2*rho-1)/(rho+1);           alphaf=rho/(rho+1);
        gamma=0.5+alphaf-alpham;            beta=0.25*(gamma+0.5)^2;
        

        % CORREGIR VELOCIDADES INICIALES, VER DE CARDONA, VA LA VELOCIDAD ANGULAR DEL PASO ANTERIOR ! PROBABLEMENTE POR ESO LA CONVERGENCIA CON AMORTIGUAMIENTO
        % ES TAN MALA
        % TAMBIEN ES IMPORTANTE PROYECTAR FUERZAS EN ALGORITMOS MULTISTEP Y NO LO ESTOY HACIENDO. REVISAR!
            A0(X.sfdof,X.tstep) = A0(X.sfdof,X.tstep) + ((1-alphaf)/(1-alpham)) * Acc(X.sfdof,X.tstep);
            
            A0(X.sfdof,X.tstep+1)=(1/(1-alpham))*(alphaf*Acc(X.sfdof,X.tstep)-alpham*A0(X.sfdof,X.tstep)); 
            
            Ui(X.sfdof,X.tstep+1) = dt*Vel(X.sfdof,X.tstep) + ((dt^2)/2)*(1-2*beta)*A0(X.sfdof,X.tstep) + (dt^2)*beta*A0(X.sfdof,X.tstep+1); 
            Ui(X.slmdof,X.tstep+1)= 0;  
            Ui(X.elidof,X.tstep+1) = FE.bc (2,:)';     %IMPOSITION OF DISPLACEMENTS IN INCREMENTAL FORM 
                        
            U(X.sfdof,X.tstep+1) = U(X.sfdof,X.tstep) + Ui(X.sfdof,X.tstep+1);
            U(X.slmdof,X.tstep+1)= 0;            
            U(X.elidof,X.tstep+1) = U(X.elidof,X.tstep) + Ui(X.elidof,X.tstep+1);     %Impose incremental displacements
            
            Vel(X.sfdof,X.tstep+1) = Vel(X.sfdof,X.tstep) + dt*(1-gamma)*A0(X.sfdof,X.tstep) + dt*gamma*A0(X.sfdof,X.tstep+1);
            %Vel(FE.ic(1,:),X.tstep+1) = FE.ic(2,:) ;  %Velocidad impuesta;
                   
            Acc(X.sfdof,X.tstep+1) = 0; 
            

    end
    
    %------------------------------------------------------------------------------
    %                                 ALPHA DE ARNOLD
    %------------------------------------------------------------------------------    
    if method ==2  
        rho = FE.step{2,2}(2);
        alpham=(2*rho-1)/(rho+1);           alphaf=rho/(rho+1);
        gamma=0.5+alphaf-alpham;            beta=0.25*(gamma+0.5)^2;
        
        Ui(X.sfdof,X.tstep+1) =  dt*Vel(X.sfdof,X.tstep) + (dt^2)*(0.5-beta)*A0(X.sfdof,X.tstep);
        U(X.slmdof,X.tstep+1)=0;

        Vel(X.sfdof,X.tstep+1) = Vel(X.sfdof,X.tstep) + dt*(1-gamma)*A0(X.sfdof,X.tstep);

        A0(X.sfdof,X.tstep+1)=(1/(1-alpham))*(alphaf*Acc(X.sfdof,X.tstep)-alpham*A0(X.sfdof,X.tstep));

        Ui(X.sfdof,X.tstep+1) = U(X.sfdof,X.tstep+1) + (dt^2)*beta*A0(X.sfdof,X.tstep+1);

        Vel(X.sfdof,X.tstep+1) = Vel(X.sfdof,X.tstep+1) + dt*gamma*A0(X.sfdof,X.tstep+1);

        Acc(X.sfdof,X.tstep+1) =  0; % ojoooo, cero aceleracion inicial no indice que Acc(tstep+1)=0, sino que el incremento es cero.

        %%OJOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        A0(X.sfdof,X.tstep+1)=(1/(1-alpham))*(alphaf*Acc(X.sfdof,X.tstep)-alpham*A0(X.sfdof,X.tstep)); 
    end



    %------------------------------------------------------------------------------
    %                                 NEWMARK UL
    %---------------------------------------------------------------------
    if method == 4  % PREDICTOR NEWMARK UL METHOD 
             
        Ui(X.sfdof,X.tstep+1) =  dt*Vel(X.sfdof,X.tstep) + ((dt^2)/4)*Acc(X.sfdof,X.tstep);
        U(X.sfdof,X.tstep+1) = U(X.sfdof,X.tstep) + Ui(X.sfdof,X.tstep+1);
        % Imposition of Geometrical Boundary Conditions
        Ui(X.elidof,X.tstep+1) = FE.bc (2,:)';     %Impose incremental displacements     
        U(X.elidof,X.tstep+1) = U(X.elidof,X.tstep) + Ui(X.elidof,X.tstep+1);     %Impose incremental displacements
        
        Vel(X.sfdof,X.tstep+1) = Vel(X.sfdof,X.tstep) + dt*0.5*Acc(X.sfdof,X.tstep);
        
        Acc(X.sfdof,X.tstep+1) = 0;
                
    end

% ------- DYNAMIC TIME EVOLUTION INFORMATION ----------------
    [RT,Ome,Gam] = rotupdd(X,RT,Ui,Vel,Acc,Ome,Gam,'pred');
    
    X.T(1,X.tstep+1) = X.tstep; 
    X.T(2,X.tstep+1) = X.tstep*FE.time(3);
    X.T(3,X.tstep+1) = FE.time(3);

end



    