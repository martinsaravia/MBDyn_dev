% --------------------------------------------------------
%            TRASLATIONAL MASS ELEMENT
% --------------------------------------------------------


function [X,gdof,m,kt,kg,fie,fme, fde,fae] = masstr(FE,X,E,el,Acc)


%----------------------------------------------------------------------------------------------- 
nn1 = E(el,4);  % node where mass is located  
mag = E(el,5);  % Mass Magnitude
edof = 3;
%-----------------------------------------------------------------------------------------------
gdof = X.CTable(nn1,1:edof);
%-----------------------------------------------------------------------------------------------

kt = zeros (edof,edof);  kg=kt; 
fie = zeros(edof,1); fde=fie; fae=fie; 


m = mag*eye(3,3);

% -------------------------------- DISTRIBUTED ELEMENT FORCES----------------------------------------
if isfield(FE,'bforce')==1
    fde = mag * FE.bforce(1:3)';  % Gravity Force
end

% -------------------------------- INCREMENTAL INERTIA FORCES----------------------------------------
fme=zeros(3,1);
if strcmp(FE.step{1,3},'DYNAMIC')==1
    DAccel    = Acc(gdof,X.tstep+1) - Acc(gdof,X.tstep);       
    fme = m * DAccel ; % Linear part of INCREMENTAL inertia forces with respect to accelerations
end