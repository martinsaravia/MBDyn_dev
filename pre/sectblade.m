%%==========================================================================
%                        BLADE SECTION PROPERTIES
%
% Computation of blade properties according to blade chord length
% Note: Modifies D matrix and assign and redifines E(:,3) to enter the
% location of element sectional propeties on D
%==========================================================================

function [X,D,DM,wf,E,count]=sectblade(X,N,FE,D,DM,E,count,j)
  

%-------------------------------------------------------------------------------------------------------------------------------------     
     nidx=j*2-1;                      didx=j*2;                    
     secfile = strcat(FE.sec{nidx,2},'.sinp');   %secfile = section name  
%-------------------------------------------------------------------------------------------------------------------------------------      
     setname = FE.sec{nidx,4};
     setidx = strmatch(setname,X.elset(:,1),'exact');  
     setels = [X.elset{setidx,2}];
     
     if isfield(FE,'rotor')==1
     rootnode = FE.rotor{didx,1}(1);     % Read root node from rotor card
     else
     disp ('Blade Section must be acoompanied by ROTOR card in order to define the rootnode')   
     end
     xroot = N(rootnode,2:4);   
        
     bladedata = FE.sec{didx,1};
%-------------------------------------------------------------------------------------------------------------------------------------          


if isempty(setels)==0
    
    for ii=1:length(setels)
        el=setels(ii); %elements that form the blades 
        nn1=E(el,4);                 nn2=E(el,5); 
        xn1 = N(nn1,2:4);            xn2 = N(nn2,2:4);
        xc = 0.5*xn1 + 0.5*xn2;                      % Element Center Position
        rpos = norm(xc - xroot);    % Position relative to root
        chord = interp1( bladedata(:,1) , bladedata(:,3) , rpos, 'cubic' );   % Spline Interpolated Chord 
        
        [NS,ES,D(:,:,count),DM(:,:,count),wf,info] = sectgen(secfile, chord, chord, 'outer');

                            % section index       name       
        X.Dinfo(count,1:2) = { FE.sec{nidx,1}  FE.sec{nidx,4}  };
        X.Dinfo(count,3:5) = info;
        
        E(el,3) = count;  % Update Element Section ID
        count = count+1;
    end
        
end


  



