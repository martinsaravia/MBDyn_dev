%============================================================
%                    PREPROCESSING
%
%  Martin Saravia - 13/08/2011 
%  Note: Connectivity Table Added
%  Note: Output file added (subfunction "infofile")
%  23/03/2011: Static Analysis Capabilities
%  13/06/2011: Joints and updated lagrangian elements
%  13/09/2011: Read from inps file
%  01/11/2012: Include Linear Multibody Elements
%  10/04/2013: Include angular velocites, no number to pregen
%============================================================

function [FE,X,N,E,D,DM,WF,F,FD,FI,FM,FA,SE,EE,U,Ui,DU,UT,Vel,Acc,Ome,Gam,A0,AV,RT]=pregen(FE)

X.flag = zeros (20);
X.wflag =zeros(20);
X.outputfile=strcat(FE.job,'.txt');
X.stop=0;  X.lambda=0; X.dlambda=0; X.iter=0; X.tstep=1;  


FE.bctype=1;     % CONDICION DE BORDE POR ELIMINACION
FE.ic=[1; 0 ; 0];  % CONDICION INICIAL
%% ==========================================================================
%                           JOB CARD READING
%==========================================================================
[X,FE,N,E] = readinps(X,FE);

FE.time = FE.step{2,1};
X.totsteps = (FE.time(2)-FE.time(1))/FE.time(3); 


FE.param(5) = 2 ; % Obligo a iterar x veces
FE.param(4) = 20; % Cantidad Maximas de Iteraaciones

% LONGITUD DE ARCO CORREGIR PARA FE.step (està para FE.dyn viejo)
if  strcmp(FE.step{2,2}(2),'arclth')==1 
    X.ArcL=1; 
    X.arclength=FE.step{2,1}(3); % LONGITUD DE ARCO INICIAL
    disp('CHEQUEAR LECTURA DE LONGITUD DE ARCO')
    X.predsign=1;
else
    X.ArcL=0;
end


infofile(X,FE,1)        % Information output file


%% ========================================================================
%                    TABLA DE CONECTIVIDAD 
%==========================================================================
actdof=1;  count=1;  lmscount=1; X.joints=0;
X.CTable=zeros(max(N(:,1)),12);
X.slmdof=[]; 
        

for i=1:X.els      
       
      node1=E(i,4);       node2=E(i,5);  node3=E(i,7);    elid=E(i,2);
   
      % Elementos  TL EU y UL y EUA
    if elid<=6 || elid==10;
            % Node 1 Processing
            if X.CTable(node1,1)==0;           
                X.CTable(node1,1:6)=actdof:actdof+5;    actdof=actdof+6;  count=count+1;
            end       
            % Node 2 Processing
            if X.CTable(node2,1)==0;           
                X.CTable(node2,1:6)=actdof:actdof+5;   actdof=actdof+6;   count=count+1;
            end
    end
    
      % ELEMENTOS LINEALES MULTICUERPO
     if elid==8 || elid==9;
            % Node 1 Processing
            if X.CTable(node1,1)==0;           
                X.CTable(node1,1:6)=actdof:actdof+5;    actdof=actdof+6;  count=count+1;
            end       
            % Node 2 Processing
            if X.CTable(node2,1)==0;           
                X.CTable(node2,1:6)=actdof:actdof+5;   actdof=actdof+6;   count=count+1;
            end
             % Node 3 Processing
            if X.CTable(node3,1)==0;           
                X.CTable(node3,1:6)=actdof:actdof+5;   actdof=actdof+6;   count=count+1;
            end
    end
    
    
    
    %   ELEMENTOS JOINT
    if elid>=80 && elid<=99; 
        X.joints=X.joints+1;
        if X.CTable(node1,1)==0;    % Node 1 Processing        
            X.CTable(node1,1:6)=actdof:actdof+5;    actdof=actdof+6;  count=count+1;
        end       
        if X.CTable(node2,1)==0;   % Node 2 Processing         
           X.CTable(node2,1:6)=actdof:actdof+5;    actdof=actdof+6;  count=count+1;
        end
        
        % SPHERICAL
         if elid==90 || elid==80
            if X.CTable(node2,7)==0; % Lagrange multipliers
              X.CTable(node2,7:9)=actdof:actdof+2;   X.slmdof=[X.slmdof actdof:actdof+2];
              actdof=actdof+3;  lmscount=lmscount+3; % Lagrange Multipliers
            end
         end
         % REVOLUTE  y CYLINDRICAL
         if elid==91 || elid==81 
            if X.CTable(node2,7)==0; % Lagrange multipliers
              X.CTable(node2,7:11)=actdof:actdof+4;   X.slmdof=[X.slmdof actdof:actdof+4];
              actdof=actdof+5;  lmscount=lmscount+5; % Lagrange Multipliers
            end
         end
    end
end

% Inactive nodes check
X.inactnodes=0;     X.rotdofs = [];
count1=0; count2=0;
 for i=1:X.nds
     if X.CTable(i,1)==0 % Si no hay GDL asignado
        count1=count1+1;
        X.inactnodes(1,count1)=i;
     else
        count2=count2+1;
        X.actnodes(1,count2)=i;
        X.rotdofs = [ X.rotdofs X.CTable(i,4:6) ];
     end
 end
 

 X.ands =length(X.actnodes);
 
 
 
 
 
 
 

%% ==========================================================================
%                        NODAL  REFERENCE TRIADS and Rotation Matrices
%==========================================================================
EE = zeros(9,X.els,X.totsteps+1);  % Strain Vector
SE = zeros(9,X.els,X.totsteps+1);  % StessVector
AV = zeros(X.els*3,6,X.totsteps+1);
RT = zeros(X.nds*3,3,4);  RT(1:3:X.nds*3,1,:) = 1;  RT(2:3:X.nds*3,2,:) = 1;  RT(3:3:X.nds*3,3,:) = 1;   % Initial Rotation (all init to ones)

for el=1:X.els
    
    elid = E(el,2);      loc=el*3-2:el*3;
       
% ------ TOTAL AND UPDATED LAGRANGIAN ELEMENTS ----------------    
    if elid==1 || elid==3  || elid==4 || elid==5 || elid==6 || elid==8 ||  elid==9 ||  elid==10;    % Elemento viga TL y UL OJO que los lineales tambien estan
        
       if elid==3  || elid==5   || elid==9 || elid==10     % BEAM ELEMENTS
           AV(loc,1:3,1)=triads(N,E,el);       AV(loc,4:6,1)=AV(loc,1:3,1);  % Reference Triad at time step 1
       end
       
      if elid==4 || elid==6 || elid==8      % BLADE ELEMENTS 
        AV(loc,1:3,1)= triadsb(FE,N,E,el);  AV(loc,4:6,1)=AV(loc,1:3,1);   % Reference Triad considering pitch angle
      end
       
       
       % ojo que para los lineales calculo las tensiones iniciales igual que para los no lineales
       nn1 = E(el,4);         nn2 = E(el,5);   
       longel=sqrt( (N(nn2,2)-N(nn1,2))^2 + (N(nn2,3)-N(nn1,3))^2 + (N(nn2,4)-N(nn1,4))^2);
       N1=0.5;         N2=0.5;           dN1=-1/longel;    dN2=1/longel;
 
       X0n1(1:3,1) = N(nn1,2:4)' ;        X0n2(1:3,1) = N(nn2,2:4)' ;             % Posicion inicial
       An1 = AV(el*3-2:el*3,1:3,1);          An2 = AV(el*3-2:el*3,4:6,1);  % Undeformed Configuration Triad (X.A is Reference Configuratio Triad)
       dX0 = dN1*X0n1 + dN2*X0n2;
       Ai = N1*An1 + N2*An2;        Ai2=Ai(:,2);       Ai3=Ai(:,3);
       dA = dN1*An1 + dN2*An2;      dA2=dA(:,2);       dA3=dA(:,3);
 
       EE(1,el,1) = dX0'*dX0;   EE(2,el,1) = dX0'*dA3;  EE(3,el,1) = dX0'*dA2;    EE(4,el,1) = dX0'*Ai2;    EE(5,el,1) = dX0'*Ai3; 
       EE(6,el,1) = dA2'*Ai3;   EE(7,el,1) = dA2'*dA2;  EE(8,el,1) = dA3'*dA3;    EE(9,el,1) = dA2'*dA3;
       
       if elid==8 || elid==9 || elid==10
        EE(:,el,1)=0;
       end
       
    end
    
% ------ EULERIAN ELEMENTS ----------------    
    if elid==2; % Elemento viga EU
       A(loc,1:3)=Ainit(N,E,el);  %#ok<AGROW>
       X.e(loc,1:3)= A(loc,1:3); % Terna en punto de integracion para elemento euleriano
       X.de(loc,1:3)=zeros(3,3);
       X.dA(loc,1:3)=zeros(3,3); 
       disp('corregir ternas y demas yerbas del elemento EU (pregen.m)')
    end
    
    if elid>=80 && elid<=89; % JOINTS TL
       AV(loc,1:3,1)=dia3(1);    AV(loc,4:6,1)=dia3(1);
   end

    if elid>=90 && elid<=99; % JOINTS UL
%      AV(loc,1:3,1)=dia3(1);    AV(loc,4:6,1)=dia3(1); %% Forma global
       AV(loc,1:3,1)=triads(N,E,el);       AV(loc,4:6,1)=AV(loc,1:3,1);  % Reference Triad at time step 1        
    end
   
end

X.sdof=actdof-1;                         % Total Degrees of Freedom of the Model (LMs indluded)
X.sfdof=1:X.sdof;
X.sfdof(X.slmdof)=[]; % Fisical DOF (no Lagrange Multipliers)
X.slms=lmscount-1;                       % Total Lagrange Multipliers of the systems


%% ==========================================================================
%                      INICIALIZACIONES 
%==========================================================================     

%% -----------------------BASE DE TIEMPO-------------------------------
X.conv=0;  % Convergence Flag
X.T = zeros(3,1);                 % Matriz de Tiempo

%% -----------------------DESPLAZAMIENTOS-------------------------------
U = zeros (X.sdof,X.totsteps+1);           % Desplazamientos (pongo 2 columnas para calcular KT Lineal 
Ui = zeros (X.sdof,X.totsteps+1);    %Desplazamiento Incrementales
UT = zeros(X.sdof,1,X.totsteps+1);
DU = zeros (X.sdof,X.totsteps+1);

if strcmp(FE.step{1,3},'DYNAMIC')==1 
    Vel = zeros (X.sdof,X.totsteps+1);           % Desplazamientos
    Acc = zeros (X.sdof,X.totsteps+1);  
    Ome = zeros (X.sdof,X.totsteps+1);           % Desplazamientos
    Gam = zeros (X.sdof,X.totsteps+1);      
    
    A0  = Acc;
    Vel (FE.ic(1,:),1) = FE.ic(2,:);
else
    Vel = 0;           % Desplazamientos
    Acc = 0;  
    A0  = 0;
    Ome = 0;
    Gam = 0;
end


%% ==========================================================================
%                      CONDICIONES DE BORDE  
%==========================================================================
idxdof=0;
bconds = length(FE.bcond(:,1));
for i=1:bconds
    bcdata = FE.bcond{i,3};
    for j=1: (length(FE.bcond{i,3}(:,1))) % Process data lines for bcond i
        node = bcdata(j,1);
        rldofs = bcdata(j,2) : bcdata(j,3); % restricted local dofs
        bcval = bcdata(j,4);
        for k=1:length(rldofs)
            idxdof=idxdof+1;
            FE.bc(1,idxdof) = X.CTable(node,rldofs(k));
            FE.bc(2,idxdof) = bcval;
        end
    end
end
     
X.elidof=FE.bc(1,:);
X.rdof=length(FE.bc(1,:));
X.srdof=X.sdof-X.rdof;
X.adof=1:X.sdof; X.adof(X.elidof)=[];

% Condiciones de borde Naturales
F = zeros(X.sdof,X.totsteps+1); 
FD = zeros(X.sdof,X.totsteps+1);
FI = zeros(X.sdof,X.totsteps+1); 
FM = zeros(X.sdof,X.totsteps+1); 
FA = zeros(X.sdof,X.totsteps+1);

%% ==========================================================================
%                            MATRICES CONSTITUTIVAS
%==========================================================================
D = zeros(9,9,X.secs);     DM = zeros(6,6,X.secs);  % Estimated PreAlloctation 
count=1;  % Counter for dimension of D and DM

for j=1:X.secs  
    
    if isempty(FE.sec{j*2-1,2})==0  
    
        if strcmp(FE.sec{j*2-1,3},'BEAM')==1;       [X,D,DM,WF,E,count] = sectbeam(X,FE,D,DM,E,count,j);        end

        if strcmp(FE.sec{j*2-1,3},'BLADE')==1;      [X,D,DM,WF,E,count] = sectblade(X,N,FE,D,DM,E,count,j);     end
        
        if strcmp(FE.sec{j*2-1,3},'GRAL')==1;       [X,D,DM,WF,E,count] = sectgral(X,N,FE,D,DM,E,count,j);      end
        
    end

end
% Write info to outputfile.
% infofile(X,FE,2)


if X.flag(10) == 1;
    X.aero = zeros(X.els,10,X.totsteps+1); 
end









function infofile(X,FE,n)

    switch (n)

    %  SALIDA A ARCHIVO

        case 1
        file_1 = fopen(X.outputfile,'w');
        % Write
        fprintf(file_1,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n');
        fprintf(file_1,'x                                                      x\n');
        fprintf(file_1,'x                MULTIBODY SBEAM                       x\n');
        fprintf(file_1,'x                                                      x\n');
        fprintf(file_1,'x              C.M. Saravia - 2011                     x\n');
        fprintf(file_1,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n\n');
        fprintf(file_1,'Job Name: %s - ',FE.job);
        fprintf(file_1,'Instant: %s  \n\n\n',datestr(now));


        case 2
        file_1 = fopen(X.outputfile,'a');    
        fprintf(file_1,'MODEL INFORMATION      \n\n');
        fprintf(file_1,' Total Elements in the Model      = %d\n', X.els);
        fprintf(file_1,' Total Nodes in the Model         = %d\n', X.nds);
        fprintf(file_1,' Total DOF in the Model           = %d\n', X.sdof);
        fprintf(file_1,' Total Joints in the Model        = %d\n', X.joints);
        fprintf(file_1,' Total Lagrange Multipliers       = %d\n\n', X.slms);

        fprintf(file_1,' MATERIAL PROPERTIES \n\n');
        fprintf(file_1,' Constitutive Law                 = %s \n', FE.mat{1});   
        if strcmp(FE.mat{1},'isot')==1
        fprintf(file_1,' Density                          = %.1f \n', FE.mat{2});        
        fprintf(file_1,' Elasticity Modulus E11           = %.1f \n', FE.mat{3});
        fprintf(file_1,' Poisson Ratio                    = %.1f \n', FE.mat{4}); 
        end
        if strcmp(FE.mat{1},'comp')==1
        fprintf(file_1,' Density                          = %.1f \n', FE.mat{2});        
        fprintf(file_1,' Elasticity Modulus E11           = %.1f \n', FE.mat{3});
        fprintf(file_1,' Elasticity Modulus E22           = %.1f \n', FE.mat{4});
        fprintf(file_1,' Shear Modulus G12                = %.1f \n', FE.mat{5});
        fprintf(file_1,' Shear Modulus G23                = %.1f \n', FE.mat{6});
        fprintf(file_1,' Poisson Ratio                    = %.2f \n\n', FE.mat{9});
        end 

        fprintf(file_1,' SECTIONAL PROPERTIES \n\n');
        fprintf(file_1,' Section Type                     = %s \n', FE.sec{1});
        fprintf(file_1,' Section Name                     = %s \n', FE.sec{2});  
        fprintf(file_1,' Thickness                        = %.1f \n\n\n', FE.sec{4});

         if strcmp(FE.step{1,3},'DYNAMIC')==1
        fprintf(file_1,' TIME INTEGRATION INFO \n\n');
        fprintf(file_1,' Total Time                       = %g \n', FE.time(2));
        fprintf(file_1,' Initial Time Step                = %g \n\n', FE.time(3));
         end

        if FE.step{2,2}(1)==2 
            fprintf(file_1,'WARNING !!!  Using the Newmark integrator for a Multibody System \n\n\n');
        end


        % Check for inactive nodes.
        if X.inactnodes(1)>=1
        fprintf(file_1,'WARNING !!! Inactive nodes:       = ');
        fprintf(file_1,'%d ', X.inactnodes);
        fprintf(file_1,'\n\n');
        %     error('xxxxx  There are inactive nodes in the model !!!  xxxxx')
        end


        fprintf(file_1,'Preprocessing Completed in... %.2f sec. \n\n',toc);

    end

fclose(file_1);
