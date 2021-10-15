%============================================================
%  PURPOSE: Compute tangent stiffness matrix
%  IN: diretors,displacements, etc
%  OUT: tangent, new directors, stresses,internal forces
%
%  Martin Saravia - 07-04-2010
% 26/10/2010 -> Agrego distintos tipos de elementos, joints y demas
%============================================================

function [X,M,K,KG,F,FD,FI,FM,FA,SE,EE,AV]=tanglin(FE,X,E,N,U,Ui,Vel,Acc,DU,D,DM,AV,F,FD,FI,FM,FA,SE,EE,RT)


M=zeros(X.sdof,X.sdof);
K=zeros(X.sdof,X.sdof);
KG=zeros(X.sdof,X.sdof);
FI(:,X.tstep+1)=0; % Mando a cero la fuerza interna total para ensamblar
FMinc=zeros(X.sdof,1); % Mando a cero la fuerza de inercia total para ensamblar
FM(:,X.tstep+1)=0; % Mando a cero la fuerza de inercia total para ensamblar
FD(:,X.tstep+1)=0; % Mando a cero la fuerza distribuida total para ensamblar
FA(:,X.tstep+1)=0; % Mando a cero la fuerza aerodinamica total para ensamblar

for el=1:X.els

    elid=E(el,2);
    
    if elid==1   % ELEMENTO SBEAM  TOTAL LAGRANGIANO
        [X,gdof,m,kt,kg,fie,fde,fme,SE,EE,AV] = sbeamtl(FE,X,E,N,D,DM,U,Ui,Vel,Acc,FD,FI,FM,SE,EE,el,AV);
    end
    
    if elid==2   % ELEMENTO SBEAM  EULERIANO      
        [m,kt,kg,fie,gdof,X] = sbeameu (X,E,N,U,DU,A,el,D,DM);
    end
    
    if elid==3  % ELEMENTO SBEAM  LAGRANGIANO ACTUALIZADO     
        [X,gdof,m,kt,kg,fie,fde,fme,SE,EE,AV] = sbeamul(FE,X,E,N,D,DM,U,Ui,Vel,Acc,FD,FI,FM,SE,EE,el,AV,RT);
    end
    
    if  elid==4  % ELEMENTO BLADE UL    
        [X,gdof,m,kt,kg,fie,fde,fme,fae,SE,EE,AV] = scp2ul(FE,X,E,N,D,DM,U,Ui,Vel,Acc,FD,FI,FM,SE,EE,el,AV,RT);
    end
    
    if elid>=80 && elid<=89; % JOINTS TL
        slavenode=E(el,5);
        lmdof=X.CTable(slavenode,7:9);
        [m,kt,kg,fie,fde,fme,X,AV] = jointtl (X,E,AV,elid,el,U,lmdof); 
        disp('arreglar fme y fde y AV proque no estan definidos para este joint tl')
        gdof=[ X.CTable(nn1,1:6) X.CTable(nn2,1:9)];
    end
    
    if elid==90     % ELEMENTO SPHERICAL JOINT UL
        slavenode=E(el,5);
        lmdof=X.CTable(slavenode,7:9);
        [m,kt,kg,fie,fde,fme,gdof,X,AV,RT] = jointul (X,E,elid,el,Ui,lmdof,AV,RT);
        fde=zeros(12,1); %No distributed force
    end

    if  elid==99 || elid==91      % ELEMENTO REVOLUTE JOINT UL
        slavenode=E(el,5);
        lmdof=X.CTable(slavenode,7:11);
        [m,kt,kg,fie,fde,fme,gdof,X,AV,RT] = jointul (X,E,elid,el,Ui,lmdof,AV,RT); 
    end
    
    
    %============================================================
    %                         ENSAMBLE  
    %============================================================
    for i=1:length(gdof)
      ii=gdof(i);
      FI(ii,X.tstep+1)  = fie(i,1) + FI(ii,X.tstep+1); % Ensamble de la fuerza interna incrementalt
      FMinc(ii,1)  = fme(i,1) + FMinc(ii,1); % Ensamble de la fuerzas de inercia del elemento (NO NODALES)
%       FM(ii,X.tstep+1)  = fme(i,1) + FM(ii,X.tstep+1); % Ensamble de la fuerzas de inercia del elemento (NO NODALES)
      FD(ii,X.tstep+1)  = fde(i,1) + FD(ii,X.tstep+1); % Ensamble de la fuerzas externas del elemento (NO NODALES)
      if elid==4; % Add aerodynamic force
        FA(ii,X.tstep+1)  = fae(i,1) + FA(ii,X.tstep+1); % Ensamble de la fuerzas externas del elemento (NO NODALES)
      end
      
      if X.tstep==1
          for j=1:length(gdof)     % Loop de gld j-esimo local
            jj=gdof(j);           % GLD global correspondiente al local j-esimo
            M(ii,jj) = m(i,j) + M(ii,jj);
            K(ii,jj) = kt(i,j) + K(ii,jj);
          end
      end
    end  
    
   
    
end

FM(:,X.tstep+1) = FM(:,X.tstep) + FMinc(:,1);  % INCREMENTAL UPDATE OF INERTIA FORCES

F = extforce (FE,X,F,U,Ui,RT);   % Update of External Forces

if X.tstep ==1
    [M,K,KG,F,FI] = bcond(FE,X,F,FI,M,K,KG,Ui);
end


