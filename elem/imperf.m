%============================================================
%  PURPOSE: Add imperfection to the system
%  IN: start and stop flags, Global displacements vector, mode
%  OUT: Imperfection
%
%  Martin Saravia - 20-04-2010(1.5.1) 
%  Modified : 23-03-2011
% 
%============================================================
function [X,U,DU,N,A]=imperf(X,M,K,U,DU,N,E,st1,st2,type,modes,amp)

[X.impvec,X.impfrq]=vmodes(M,K,modes);

% Imperfection in Displacement vector
if strcmp(type,'disp')==1 
  if X.tstep>=st1 && X.tstep<=st2
    U (X.adof,X.tstep+1) = U (X.adof,X.tstep+1) + amp * X.impvec (:,modes);
    DU (X.adof,1) = DU (X.adof,1) + amp * X.impvec (:,modes);
  end
end



% Imperfection in geometry
if strcmp(type,'geom')==1 
    imperf=zeros(X.sdof,1);
    imperf(X.adof,1) = X.impvec;
    N (:,2) = N (:,2) + amp * imperf (1:6:X.sdof,1);
    N (:,3) = N (:,3) + amp * imperf (2:6:X.sdof,1);
    N (:,4) = N (:,4) + amp * imperf (3:6:X.sdof,1);
    
    %% ==========================================================================
    %                 MODIFICAMOS LA TERNAS NODALES 
    %==========================================================================
    % Ternas de referencia nodales
    for el=1:X.els
        if E(el,4)==1; % Elemento viga
           A(el*3-2:el*3,1:3)=Ainit(N,E,el); 
        end
    end
    A(:,4:6)=A(:,1:3);
    for el=1:X.els
        if E(el,4)==99; % Elemento joint
           A(el*3-2:el*3,1:3)=A((el-1)*3-2:(el-1)*3,4:6);
           A(el*3-2:el*3,4:6)=A((el+1)*3-2:(el+1)*3,1:3);
        end
    end

end




    

    