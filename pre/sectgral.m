%============================================================
%                GEN SECTION PROPERTIES
%
% Reading an input .stn file the routine obtains its sectional
% stiffness matrix acording to Saravia-2012
%============================================================
function [X,D,DM,wf,E,count] = sectgral(X,N,FE,D,DM,E,count,j)

%------------------------------------------------------------------------------
nidx=j*2-1;                      didx=j*2;    

secfile = strcat(FE.sec{nidx,2},'.sinp');   %secfile = section name  

% secfile = [ FE.rundir  FE.sec{nidx,2} '.sinp'];   %secfile = section name 

% formu = 'NOLINEAR'; %
formu = FE.sec{nidx,7};

%-----------------------------------------------------------------------------   

setname = FE.sec{nidx,4};
setidx = strmatch(setname,X.elset(:,1),'exact');  
setels = [X.elset{setidx,2}];

%-------------------------------------------------------------------------------------------------------------------------------------
if strcmp(formu,'LINEAR')
    [NS,ES,D(1:6,1:6,count),DM(:,:,count),LS,MS,wf,info] = sectgen(secfile, 1, 1, 'middle','NO',formu);
     if X.wflag(3)==0
         disp ('INFORMATION! Using 6x6 Sectional Stiffness')
         X.wflag(3)=1;
     end

else
    [NS,ES,D(1:9,1:9,count),DM(:,:,count),LS,MS,wf,info] = sectgen(secfile, 1, 1, 'middle','NO',formu);
     if X.wflag(3)==0
         disp ('INFORMATION! Using 9x9 Sectional Stiffness')
         X.wflag(3)=1;
     end
end

X.Dinfo(count,1:2) = { FE.sec{nidx,1}  FE.sec{nidx,4}  };
X.Dinfo(count,3:5) = info;
E(setels,3) = count;  % Update Element Section location
E(setels,8) = j;  % Update Element Section ID
count = count+1;

%---------------------------------------------------------------
%                        WARNING FLAGS
%---------------------------------------------------------------        
 if X.wflag(2)==0
     disp ('WARNING! Correct sectgen.m Mass Centroid calculation')
     X.wflag(2)=1;
 end  

        
        
        
        
        
        
        
        
        