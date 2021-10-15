
function [wf,dwfs] = warpfcn (NS,ES,LS,MS,elset) 

els = length(ES(:,1));
nds = length(NS(:,1));
sets = length (elset(:,1));

%-----------------------------------------------
%            DEGREES OF FREEDOM TABLE
%-----------------------------------------------
DOF = zeros (nds,2);
DOF(:,1) = NS(:,1)' ;
gdof=1;
for j=1:sets
    setels = elset{j,2} ;  
    for i=1:length (setels)
       el=setels(i);    
        for ind=1:2
            nodo=ES(el,ind+1);
            if DOF(nodo,2)==0; 
                DOF(nodo,2)=gdof;               
                gdof=gdof+1;
            end
        end

    end
end
dofs = gdof-1;
%-----------------------------------------------
%           STIFFNESS AND LOAD VECTOR
%-----------------------------------------------
K = zeros(dofs,dofs);
F=zeros(dofs,1);
W=zeros(dofs,1);

for j=1:sets
    setels = elset{j,2} ;    
    setlamid = strmatch( elset{j,3}, LS(:,1), 'exact');
    
    for i=1:length (setels)
      
        el=setels(i);    n1=ES(el,2);   n2=ES(el,3);
        edof = [DOF(n1,2) DOF(n2,2)];
        p1=NS(n1,3:4);  p2=NS(n2,3:4);   
        Le = norm(p2-p1);
        
        y=NS(n1,3); z=NS(n1,4);     % La coordenada rn de cualuier nodo es la misma porque el segmento es recto
        dys=(NS(n2,3)-NS(n1,3))/Le; % ds=LE;
        dzs=(NS(n2,4)-NS(n1,4))/Le;
        rn = y*dzs-z*dys;
        
        [CL,wseg,thseg] = lamseg(LS,MS,setlamid,'BRB');
        Gseg = CL(2,2)/thseg; % G * th
   
        kt = (Gseg*thseg/Le)* [ 1  -1  ;
                               -1   1  ];

        fe = (Gseg*(rn*thseg+thseg^2/2)/1) * [1 ; -1];

        if rn<=0
           disp(['WARNING: Negative normal coordinate for segment  ' num2str(el) ' during warpfcn'])   % rn debe ser siempre positivo (sino esta mal la orientacion)

        end

        for lf=1:2  % Local File Loop
            gf=edof(lf);
            F(gf,1)=fe(lf,1) + F(gf,1);
            for lc=1:2  % Local Column Loop
                gc=edof(lc);
                K(gf,gc) = kt(lf,lc) + K(gf,gc);
            end
        end
    end
end

% save('K','K')

% [V,D] = eig(K(2:dofs,2:dofs))

% [K,F]=penal(K,F,[1 0]);
% wf=K\F;
% wf=wf-sum(wf)/length(wf); % Total zero displacement
 
% --------------------------------------------------
%   SOLVE THE PROBLEM BY ELIMINATION OF CONSTRAINTS
% --------------------------------------------------

wf(2:dofs,1)=K(2:dofs,2:dofs)\F(2:dofs,1);

  
wf = wf - sum(wf)/length(wf); % Impose constraint of total zero displacement

% wf
% F   
% JT = F' * wf



% --------------------------------------------------
%  CALCULATE DERIVATIVE OF W RESPECT TO s
% --------------------------------------------------

dwfs = zeros(els,2);
for j=1:sets
    setels = elset{j,2} ;      
    for i=1:length (setels)     
        el=setels(i);    n1=ES(el,2);   n2=ES(el,3);
        p1=NS(n1,3:4);   p2=NS(n2,3:4);   
        Le = norm(p2-p1);
        dwfs(el,1) = el;
        dwfs(el,2) = (wf(DOF(n2,2))-wf(DOF(n1,2))) / Le;
    end
end




