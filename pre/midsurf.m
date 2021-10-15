%%===========================================================
%                 Extraction of midsurface
%
% In: CLOCKWISE ORDERED outer surface (matrix NS), average thickness, direction
% Out: midsurface.
%===========================================================


 function  NSm = midsurf(NS,ES,LS,MS,surf)

 
 nds =length(NS(:,1));
 NSm = zeros (nds,4); % Midsurface Coordinates
 NSm(:,1) = NS(:,1);  % Node numbering
 

 %  Choose Inner or outer direction
 if strcmp(surf,'outer')==1
    iof = -1; %Inner or outer direction
 elseif strcmp(surf,'inner')==1
    iof = 1;
 else
    iof = 0;
 end
 
 
 for i=1:nds
     
    nnodo = NS(i,1); % Node number
    
    % CARTESIAN BASIS     
    iv = [1 0 0];    jv=[0 1 0];     kv=[0 0 1]; 
     
    %ORIGINAL (NON ORTHOGONAL) BASIS OF NORMAL VECTORS
   
    [ela,col,nn] = find (ES(:,3)==nnodo); % Elemento anterior
    if isempty (ela)== 1; continue; end % Skip not used nodes
    tha = LS{ ES(ela,4), 6 }; % Espesor del elemento anterior
    
        
    [elp,col,nn] = find (ES(:,2)==nnodo); % Elemento posterior
    thp = LS{ ES(elp,4), 6 }; % Espesor de elemento posterior
    
    if nnodo == 29 || nnodo == 31 ||nnodo == 32 ||nnodo == 34 
            t = [0 tha/2 thp/2];  
    else
            t = [0 tha/2 thp/2];  % Covariant components of tv in original basis  
    end

    
%     if i==1;       
%         d2 = NS(i,2:4) - NS(nds,2:4); %By-pass problem of first node
%         d3 = NS(i+1,2:4) - NS(i,2:4);
%     elseif i==nds; 
%         d2 = NS(i,2:4) - NS(i-1,2:4); 
%         d3 = NS(1,2:4) - NS(i,2:4)  ; %By-pass problem of last node
%     else
%         d2 = NS(i,2:4) - NS(i-1,2:4);
%         d3 = NS(i+1,2:4) - NS(i,2:4);
%     end
    
    d2 = NS( ES(elp,3) , 2:4) - NS( ES(elp,2) , 2:4); % Longitud del elemento posterior vector gsub2
    d3 = NS( ES(ela,2) , 2:4) - NS( ES(ela,3) , 2:4); % Longitud del elemento posterior vector gsub3

    n2 =- [0 -d2(3)  d2(2)];        
    n3 =- [0  d3(3)  -d3(2)];  % Original Basis

    n1=[1 0 0];         n2=n2/norm(n2);         n3=n3/norm(n3);
    

    
    % RECIPROCAL BASIS
    v = n3 * f_cross(n1,n2);  % Volume or paralellepiped
    
    
    e1=cross(n2,n3)/v;  
    e2=cross(n3,n1)/v;  
    e3=cross(n1,n2)/v; % Reciprocal Basis
    
    % TRANSFORMATION
    M = [iv*e1'  jv*e1'  kv*e1';
         iv*e2'  jv*e2'  kv*e2';
         iv*e3'  jv*e3'  kv*e3'];  % Transformation matrix, Original-Cartesian

    tv = M'*t'; % Cartesian components of vector

    % NEW POSITION OF REFERENCE POINT  sum(isfinite(tv))>=3
    
    if abs(v) >= 1E-1 % Low pass value for semi aligned points
        NSm(i,2:4) = NS(nnodo,2:4) + iof * tv';
    else
         NSm(i,2:4) = NS(nnodo,2:4) + iof * (0.5*(t(2)+t(3)))*(0.5*(n2+n3)); %By-pass for colinear points (v=0) 
        disp(['WARNING: BYPASSING COLLINEAR POINTS NODE ' num2str(nnodo)]) 
    end 
     
 end

scatter (NS(:,3),NS(:,4)); hold on; axis equal
scatter (NSm(:,3),NSm(:,4),'filled')