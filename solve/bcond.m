%============================================================
%  PURPOSE: Apply BCs 
%  23-03-2011
%============================================================

function [M,C,K,KG,F,FI]=bcond(FE,X,F,FI,M,C,K,KG,U)

FE.bctype=1;

%============================================================
%  PURPOSE: Apply BCs by elimination
%  IN: Stiffness Matrix, Ext Force Vector and Int Force Vector
%  OUT: Reduced vector and matrices and Corrected Ext Force.
%
%  Martin Saravia - 27-11-2009 - 07-04-2010
%
%   NOTA: ACTUALIZAR PARA DESPLAZAMIENTOS IMPUESTOS QUE NO SEAN CERO
%============================================================

   if FE.bctype==1; 
    % Armo los vectores que contienen las columnas a eliminar de K para luego corregir por desplazamientos impuestos
%     BCK=K(X.adof,X.elidof); 
    
    % Reduction 
    M = M(X.adof,X.adof);
    C = C(X.adof,X.adof);
    K = K(X.adof,X.adof);

    if strcmp(FE.step{1},'buckle')==1
        KG = KG(X.adof,X.adof);
    end

    % External forces correction for imposed displacements. Evaluo nuevamente el vector de carga sino se actualiza incrementalmente con BCK y esta mal
    % Notar que ahora el vector de carga es variable ya que BCK es variable con la configuracion
    
    
%      F(X.adof,X.tstep+1) =  F(X.adof,X.tstep+1) - BCK * U(X.elidof,X.tstep+1); %4.43 Bathe

%% OJOOOOO--- CUANDO DEJO EL SEGUNDO TERMINO DE LA ECUACION ANTERIOR NO ANDA, SI LO SACO ANDA PERFECTO... VER PORQUE, PUEDE SR PORQUE TENIA MAL LA RUTINA
%% UPDATED LAGRANGIAN. Esta mal como lo hice, porque F se actualiza iterativamente y su sobremodifica, pero igual corrigiendolo no anda, ya lo probe.
     
   end




%============================================================
%  PURPOSE: Apply BCs by PENALIZATION
%  IN: Stiffness Matrix, Ext Force Vector and Int Force Vector
%  OUT: Reduced vector and matrices and Corrected Ext Force.
%
%  Martin Saravia - 27-11-2009 - 07-04-2010
%
%   NOTA: ACTUALIZAR PARA DESPLAZAMIENTOS IMPUESTOS QUE NO SEAN CERO
%============================================================


%ARREGLAR VETORES F Y T PARA T=tstep+1
if FE.bctype==2;
    X.Cpen=1E20;
    nrestri=length(FE.bc(1,:));
    for i=1:nrestri
        dof=FE.bc(1,i);
        K(dof,dof)=K(dof,dof)+X.Cpen;
        F(dof)=X.Cpen*FE.bc(2,i);
        if X.iter>=2
          FI(dof,X.tstep+1)=FI(dof,X.tstep+1)+X.Cpen*FE.bc(2,i);
        end
    end
end