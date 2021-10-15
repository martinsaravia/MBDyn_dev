%============================================================
%                GEN SECTION PROPERTIES
%
% Reading an input .stn file the routine obtains its sectional
% stiffness matrix acording to Saravia-2012
%============================================================
function [NS,ES,D,DM,LS,MS,wf,info] = sectgen(file, ampy, ampz, surf, cent, formu)

%% =======================================================================
%                 READ SECTION DATA FROM FILE
%=========================================================================
fid=fopen(file,'r');
node=1; elem=1; mate=1; lami=1; setid=1;
line = fgetl(fid);

while strncmp(line,'*END',4)==0 
 
        %NODES PRE
        if strfind(line,'*NODE')==1   
             line = fgetl(fid);
            while strncmp(line,'*',1)==0
                data = str2num(line);  %
                NS(data(1),:)= data;   % %SUMAMENTE IMPORTANTE, ALMACENO EN FILA IGUAL A SU NUMERO
                node = node + 1;
                line = fgetl(fid); %next line
                continue
            end
                
          % ELEMENTS PRE
        elseif strfind(line,'*ELEMENT')==1   
            
            elname=regexp(line, '(?<=TYPE=)......', 'match','once'); %#ok<NASGU> %selecciona las 4 primeras letras del nombre de elemento, once larga un string, sino es lista de strings la salida.             
            setname = regexp(line, '(?<=ELSET=)\w*', 'match','once');  
            setlami = regexp(line, '(?<=LAMINATE=)\w*', 'match','once');  
            setels=[]; %new counter for elements of the set 'setname'   
            setel=1;  %set elements counter
            
             line = fgetl(fid);
            while strncmp(line,'*',1)==0
                data = str2num(line); 
                ES(data(1),:)= data; %#ok<AGROW> %SUMAMENTE IMPORTANTE, ALMACENO EN FILA IGUAL A SU NUMERO
                setels(1,setel) = ES(data(1),1); %#ok<AGROW>
                
                setel=setel+1;
                elem = elem + 1;
                
                line = fgetl(fid); %next line
                continue
            end
            
            elset(setid,:) = { setname   setels  setlami}; %#ok<*AGROW>
            setid=setid+1;
  

                        %NODES PRE
         elseif strfind(line,'*MATERIAL')==1 
             
             matname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña

             line = fgetl(fid);   
             matlaw = regexp(line, '(?<=*)\w*', 'match','once');
             mattype = regexp(line, '(?<=TYPE=)\w*', 'match','once');
             
              line = fgetl(fid);
            while strncmp(line,'*',1)==0
                data = str2num(line); 
                MS(mate,:)= { matname  matlaw  mattype data };            
                mate = mate + 1;
                line = fgetl(fid); %next line 
                continue
            end

        
         elseif strfind(line,'*LAMINATE')==1 
             
            lamname = regexp(line, '(?<=NAME=)\w*', 'match','once');
            lammate = regexp(line, '(?<=MATERIAL=)\w*', 'match','once');
                      
            line = fgetl(fid);
            while strncmp(line,'*',1)==0
                data1 = str2num(line);  %#ok<*ST2NM> % LAMINATION SEQUENCE
                line = fgetl(fid);
                data2 = str2num(line);  % LAYER THICKNESSES
                line = fgetl(fid);
                data3 = str2num(line);  % LAYER MATERIAL: mejorar, lee numero de material, no nombre... debo conocer de antemano el nombre.
                data4 = sum(data2); % LAMINATE THICKNESS
                
                LS(lami,:)= {lamname lammate data1 data2 data3 data4}; 
                
                lami = lami + 1;
                line = fgetl(fid); %next line 
                continue
            end
        else
         line = fgetl(fid); %next line
        end
end



%% ----------------------------------------------------------------------
%             ASSIGN LAMINATE NUMBER TO ELEMENT TABLE
% ----------------------------------------------------------------------
sets = length (elset(:,1));
for j=1:sets
    setels = elset{j,2} ;       
    setlamid = strmatch( elset{j,3}, LS(:,1), 'exact'); 
    ES(setels,4) = setlamid;
end


%% ----------------------------------------------------------------------
%                CORRECTION OF NEGATIVE NORMALS - NO SE SI FUNCIONA DEL TODO BIEN
% ----------------------------------------------------------------------

sets = length (elset(:,1));
% 
for j=1:sets
    setels = elset{j,2} ;       
    for i=1:length (setels)   % Find centroid and static moments
        el=setels(i);    
        n1=ES(el,2);    n2=ES(el,3);
        p1=NS(n1,3:4);  p2=NS(n2,3:4);          
        Le = norm(p2-p1);        
        sv = (p2-p1)/Le;             % versor tangente 
        nv = [ sv(2) -sv(1) ]; 
               
        rn1 = p1*nv'; % coordenada normal del nodo 1 (de cualquiera seria lo mismo)
        rt1 = p1*sv'; % coordenada tangente del nodo 1 (de cualquiera seria lo mismo)
        
        if rn1 < 0   %%&&  rt1 < 0 CHEQUEAR ESTA
            ES(el,2)=n2;  ES(el,3)=n1;  % si el vector tangente es horario lo invierto a antihorario (2.1.4 Librescu)
            disp(['Normal for segment ' num2str(el) ' was corrected'])
        end
    end
end

        
%% ----------------------------------------------------------------------
%                 AMPLIFICATION AND OFFSET OF SURFACE
% ----------------------------------------------------------------------
%  Amplification of section coordinates
NS(:,3) = ampy * NS(:,3);
NS(:,4) = ampz * NS(:,4);

if strcmp(surf,'outer') ==1 || strcmp(surf,'inner')==1
    
    % Deletion of web elements in order to avoid their processing in the midsurface algorithm
    ESnw = ES; % Sectional Elements Table without web elements (webx must be the name of set)
    for j=1:sets
        elset(j,1)
        if strncmpi(elset(j,1),'web',3)==1
            setels = elset{j,2} ;
            ESnw(setels,2:3)=0; % Mando a cero los nodos de los elementos web (para que midsurf no los encuentre)
            disp (['INFO: Web ' elset(j,1) ' deleted from midsurface algorithm:'])
        end
    end
    
    % TRANSPORT SECTION TO MIDLINE
    disp('Reference line transported to midline...')
    NS = midsurf(NS,ESnw,LS,MS,surf); %Change Coordinates of midline
end




%% =======================================================================
%                        CONSTITUTIVE MATRIX
%=========================================================================

areat=0; Qy=0; Qz=0;   

%-----------------------------------------------------------
%                 CENTROID POSITION  
%----------------------------------------------------------- 
for j=1:sets

    setels = elset{j,2} ;       
    setlamid = strmatch( elset{j,3}, LS(:,1), 'exact');
    setlami = LS(setlamid,:);

    for i=1:length (setels)   % Find centroid and static moments

        el= setels(1,i);     
        n1=ES(el,2);           n2=ES(el,3);            
        y0=NS(n1,3);           y1=NS(n2,3);
        z0=NS(n1,4);           z1=NS(n2,4);           
        Dy=y1-y0;              Dz=z1-z0;       Ds=sqrt( (Dy)^2 + (Dz)^2 );

        thick = LS{setlamid,6};

        areai = thick*Ds;      areat = areai + areat;
        Qy = areai*(y0+y1)/2 + Qy ;   Qz = areai*(z0+z1)/2 + Qz ;  % Momentos Estaticos
    end  

end

YC = Qy/areat;    ZC = Qz/areat;

% % TRASLATION OF SECTION POLE TO THE CENTROID 
if strcmp(cent,'YES')
    NS(:,3) = NS(:,3) - YC  ;
    NS(:,4) = NS(:,4) - ZC   ;      
end

%-----------------------------------------------------------
%           CALCULATE THE WARPING FUNCTION  
%-----------------------------------------------------------  
[wf,dwfs] = warpfcn (NS,ES,LS,MS,elset); %#ok<ASGLU>


%-----------------------------------------------------------
%           MATRICES SECCIONALES D Y DM   
%-----------------------------------------------------------     


if strcmp(formu,'LINEAR') ==1
    D=zeros(6,6);  
else
    D=zeros(9,9);  
end


DM = zeros(6,6);     asect=0;     
stds=0;% MIDLINE LENGTH

for j=1:sets
    setels = elset{j,2} ;       
    setlamid = strmatch( elset{j,3}, LS(:,1), 'exact');
    
%     s0=0;   % Reinit s coordinate for every set (it should be the rigth thing to do)  
    
    for i=1:length (setels)
        s0 = 0;  %OJO, REINICIO LA COORDENADA S PARA CADA SEGMENTO
                
        el = setels(1,i);       n1=ES(el,2);      n2=ES(el,3); 
        y0=NS(n1,3);            y1=NS(n2,3);      z0=NS(n1,4);    z1=NS(n2,4);
        Dy=y1-y0     ;          Dz=z1-z0;
        Ds=sqrt( Dy^2 + Dz^2 ); s1=s0+Ds; stds=stds+Ds;
        dYs=Dy/Ds;  dZs=Dz/Ds;

        r0 = sqrt(y0^2+z0^2);   r1 = sqrt(y1^2+z1^2);  spr = 0.5*(Ds+r0+r1);
        asect = sqrt(spr*(spr-Ds)*(spr-r0)*(spr-r1)) + asect;  % Formula de Heron para area de triangulo en funcion de sus lados
        
        dwfsseg = dwfs(el,2);

        [CL,wseg,thseg] = lamseg(LS,MS,setlamid,'BRB'); %#ok<NASGU>

        if strcmp(formu,'LINEAR') ==1
            [Dseg,DMseg] = DJseg_AH_w_L (el,CL,wseg,s0,s1,Dy,Dz,Ds,y0,z0,dYs,dZs,dwfsseg);
        else
            [Dseg,DMseg] = DJseg_AH_w (el,CL,wseg,s0,s1,Dy,Dz,Ds,y0,z0,dYs,dZs,dwfsseg); 
        end
        
        D  = D  + Dseg;
        DM = DM + DMseg;
        
 
    end

    %       centroid  sectarea    massdensity
    info={ [YC ZC]      asect       DM(1,1)};
    

    
end
        
        

        
        
        
        
        
        
        