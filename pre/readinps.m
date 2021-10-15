function [X,FE,N,E] = readinps(X,FE)

filename=strcat(FE.job,FE.ext);

fid=fopen(filename,'r');

line = fgetl(fid);

E = zeros(1,6);
N = zeros(1,4);


%% ========================================================================
%                         READ NODE DATA
%==========================================================================
nodo=0; % Avanza linea para leer nodos
while strncmp(line,'*END',4)==0  
    if strncmp(line,'*NODE',4)==1       
        line = fgetl(fid);
        while strncmp(line,'*',1)==0
            nodo = nodo + 1;
            data = str2num(line); 
            N(data(1),:)= data;
            line = fgetl(fid); %next line
        end
        break
    else
        line = fgetl(fid); %adquiere linea hasta encontrar *NODE  
    end
end
X.nds=nodo;

%% ========================================================================
%                         READ ELEMENT DATA
%==========================================================================
el=0; setid=0;
while   isempty(strfind(line,'SECTION'))==1
       
     if strncmp(line,'*ELEMENT',8)==1
       
        elname=regexp(line, '(?<=TYPE=)......', 'match','once'); %selecciona las 4 primeras letras del nombre de elemento, once larga un string, sino es lista de strings la salida.             
        setname = regexp(line, '(?<=ELSET=)\w*', 'match','once');  
        setels=0; %new counter for elements of the set 'setname'   
        setel=1;  %set elements counter

        if ischar(setname)==1
           setid=setid+1;
        end
                    
       switch elname
            case 'SCB2TL';  elformu=1 ; ndsxel=3;
            case 'SCB2EU';  elformu=2 ; ndsxel=3;
            case 'SCB2UL';  elformu=3 ; ndsxel=3;
            case 'SCP2UL';  elformu=4 ; ndsxel=3;
            case 'SCB2LL';  elformu=5 ; ndsxel=3;
            case 'SCP2LL';  elformu=6 ; ndsxel=3;
            case 'SCP3LL';  elformu=8 ; ndsxel=4;  % Linear Lagrangian Multibody Element
            case 'SCB3LL';  elformu=9 ; ndsxel=4;
            case 'SSB2UL';  elformu=10; ndsxel=3;  % Small Strain UL Beam    
            case 'JOINTT';  elformu=89; ndsxel=3; 
            case 'SPHJUL';  elformu=90; ndsxel=3;
            case 'REVJUL';  elformu=91; ndsxel=3;
            case 'MASSTR';  elformu=7;  ndsxel=2;
            otherwise; disp('UNKNOWN ELEMENT TYPE');
        end
        
        line = fgetl(fid);   
        
        while strncmp(line,'*',1)==0
            el = el + 1;
            linenum=str2num(line); %#ok<*ST2NM> %linea como numero
            E(linenum(1),1)=linenum(1);                        % numero del elemento
            E(linenum(1),2)=elformu;                           % formulacion del elemento
            if elformu>=80; E(linenum(1),7)=linenum(5); end    % axis for joint elements
            
            E(linenum(1),4:(3+ndsxel))= linenum(2:(1+ndsxel)); % nodos del elemento    MEJORAR, ES UNA PORQUERIA ESTA ASIGNACION
            line = fgetl(fid); %next line
            setels(1,setel) = E(linenum(1),1); %#ok<*AGROW>
            setel=setel+1;
            
        end
           elset(setid,:) = { setname   setels}; % New element to set

        continue
     else
         line = fgetl(fid); %adquiere linea hasta encontrar *ELEMENT
     end
 
end
X.els=el;
X.elset = elset;



%% ========================================================================
%                         READ UNTIL END
%==========================================================================
% Initizalizationo of ids.
secid=0; matid=0; lamid=0; rotorid=0; windid=0; stepid=0;  bcondid=0; plotid=0; pforceid=0;
% Initizalizacion of counts
X.plotcount=0; %
while strncmp(line,'*END',4)==0 
    
        
    %================== BEAM SECTION =======================       
    if strfind(line,'*BEAM SECTION')==1 % si encuentra la palabra seccion despues del segundo caracter
        secid=secid+1; 
        %          sectype = regexp(line, '(?<=*)\w*', 'match','once'); % 4 letters of the section type
        setname = regexp(line, '(?<=ELSET=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
        secname = regexp(line, '(?<=SHAPE=)\w*', 'match','once');
        seclaw =  regexp(line, '(?<=LAW=)\w*', 'match','once');
        secmat = regexp(line, '(?<=MATERIAL=)\w*', 'match','once');
        seclam = regexp(line, '(?<=LAMINATE=)\w*', 'match','once');

        line = fgetl(fid);
        secdata = str2num(line); %read section data and transform to number

        FE.sec(2*secid-1,1:7) = { secid  secname  'BEAM'   setname  seclaw   secmat   seclam } ;
        FE.sec(2*secid,1:2) = { secdata(1) secdata(2:3)};  

    %%================== BLADE SECTION ======================= 
    elseif strfind(line,'*BLADE SECTION')==1 % si encuentra la palabra seccion despues del segundo caracter
        secid=secid+1; 
        secname = regexp(line, '(?<=NAME=)\w*', 'match','once');
        seclaw =  regexp(line, '(?<=LAW=)\w*', 'match','once');
        setname = regexp(line, '(?<=ELSET=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
        bldfile = regexp(line, '(?<=BLDFILE=)\w*', 'match','once'); % File for blade geometry

        FE.sec(2*secid-1,1:7) = { secid  secname  'BLADE'  setname  bldfile 0 seclaw} ;

        bladedata = dlmread(strcat(bldfile,'.btb')); % Read blade geometry

        if exist(strcat(secname,'.ftb'),'file')==0 % IF file doesnt exist
         foildata = [0 0 0 0 0 0]; 
         disp(['Foil Aerodynamic Table ' strcat(secname,'.ftb') ' does not exist' ])
        else
            foildata  = dlmread(strcat(secname,'.ftb')); % Aerodynamic data
        end

        FE.sec(2*secid,1:2) = { bladedata  foildata };

    %%================== GENERAL SECTION ======================= 
    elseif strfind(line,'*GRAL SECTION')==1 % si encuentra la palabra seccion despues del segundo caracter
        secid=secid+1; 
        secname = regexp(line, '(?<=NAME=)\w*', 'match','once');
        seclaw =  regexp(line, '(?<=LAW=)\w*', 'match','once');
        setname = regexp(line, '(?<=ELSET=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
        foilfile = regexp(line, '(?<=FOIL=)\w*', 'match','once');    
        bladefile = regexp(line, '(?<=BLADE=)\w*', 'match','once') ;

        if isempty (foilfile)==0,   foilfile = [foilfile '.ftb'];    foilfile = dlmread(foilfile);  end
        if isempty (bladefile)==0,  bladefile = [bladefile '.btb'];  bladefile = dlmread(bladefile); end
        FE.sec(2*secid-1,1:7) = { secid  secname  'GRAL'  setname   foilfile  bladefile seclaw}; 
        FE.sec(2*secid,1) = { 'no data lines'};

    %%================== MATERIAL LINE ======================= 
    elseif strfind(line,'*MATERIAL')==1 % si encuentra la palabra seccion despues del segundo caracter
        matid=matid+1;
        matname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña

        line = fgetl(fid);   
        matlaw = regexp(line, '(?<=*)\w*', 'match','once');
        mattype = regexp(line, '(?<=TYPE=)\w*', 'match','once');

        line = fgetl(fid);
        matdata = str2num(line); %read section data and transform to number

        FE.mat(2*matid-1,:)  = { matid   matname   matlaw  mattype  matdata};
        FE.mat(2*matid,1)  = { matdata };

    %================== LAMINATE LINE =======================       
    elseif strfind(line,'*LAMINATE')==1 % si encuentra la palabra despues del segundo caracter
         lamid=lamid+1;
         lamname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña

         line = fgetl(fid);
         lamdata = str2num(line); %read section data and transform to number

         FE.lam(2*lamid-1,:) = { lamid    lamname }; 
         FE.lam(2*lamid,1)   = { lamdata  }; 


    %================== ROTOR LINE =======================
    elseif strfind(line,'*ROTOR')==1 % si encuentra la palabra despues del segundo caracter
         rotorid=rotorid+1; if rotorid>=2 ; disp('WARNING! Only one *ROTOR is allowed'); end
         rotorname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         rotorset = regexp(line, '(?<=ELSET=)\w*', 'match','once');
         rotorstart = regexp(line, '(?<=START=)\w*', 'match','once');

         line = fgetl(fid);
                    % hubnode  hingel 
         rotorctrl1 = str2num(line); %read ROTOR data and transform to number

         if strcmp(rotorstart,'YES')==1;
             line = fgetl(fid); 
                      %      FLAG    DATA
             rotorctrl2 = [   1  str2num(line) ]; 
         else
             rotorctrl2 = 0;
         end

         FE.rotor(1,1:3) = { rotorid    rotorname   rotorset}; 
         FE.rotor(2,1:2)   = { rotorctrl1 rotorctrl2 };    

         %================== STEP LINE =======================
     elseif strfind(line,'*STEP')==1 % si encuentra la palabra despues del segundo caracter
         stepid=stepid+1;
         stepname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         steptype = regexp(line, '(?<=TYPE=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña

         % Control parameters
         stepctrl1 = regexp(line, '(?<=METHOD=)\w*', 'match','once') ;
%              regexp(line, '(?<=NUMDAMP=)\w*', 'match','once')

         stepctrl2 = str2num ( regexp(line, '(?<=NUMDAMP=)\w*\.\d*', 'match','once') );
         stepctrl3 = regexp(line, '(?<=NLGEOM=)\w*', 'match','once') ;
         stepctrl4 = str2num ( regexp(line, '(?<=STRDAMP=)\w*\.\d*', 'match','once') );
         stepctrl5 = str2num ( regexp(line, '(?<=STRFRQ=)\w*\.\d*', 'match','once') );

         if strcmp(stepctrl3,'NO') == 1;  stepctrl3 =0;  end
         if strcmp(stepctrl3,'YES') == 1;  stepctrl3 =1;  end

         if strcmp(stepctrl1,'ALPHA') == 1;  stepctrl1 =1;  end
         if strcmp(stepctrl1,'ALPHA2') == 1;  stepctrl1 =2;  end
         if strcmp(stepctrl1,'NEWMARKTL') == 1;  stepctrl1 =3;  end
         if strcmp(stepctrl1,'NEWMARKUL') == 1;  stepctrl1 =4;  end

         if isempty(stepctrl1) == 1; stepctrl1=1;  end   % DEFAULT METHOD
         if isempty(stepctrl2) == 1; stepctrl2=0.9; end % DEFAULT NUMERICAL DAMPING 
         if isempty(stepctrl3) == 1; stepctrl3=1; end % DEFAULT NLGEOM=YES
                     % METHOD   DAMPCOEF   NLGEOM
         stepctrl = [stepctrl1 stepctrl2 stepctrl3 stepctrl4 stepctrl5];

         line = fgetl(fid);
         stepdata = str2num(line); %read section data and transform to number


         FE.step(2*stepid-1,1:3) = { stepid    stepname   steptype }; 
         FE.step(2*stepid,1:2)   = { stepdata  stepctrl };

     %================== BOUNDARY LINE =======================
     elseif strfind(line,'*BOUNDARY')==1 % si encuentra la palabra despues del segundo caracter
         bcondtype = regexp(line, '(?<=TYPE=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         bcondid=bcondid+1;
             line = fgetl(fid);
             lidx = 0 ;

        while isempty(strfind(line,'*'))==1 
             lidx = lidx+1;
             bconddata(lidx,:) = str2num(line); %read section data and transform to number  
             curpos = ftell(fid); %get current position (to go back after the end of loop)
             line = fgetl(fid);
        end

        curpos = fseek(fid,curpos,'bof'); % move lo last file in order to avoid skeepint commands
        FE.bcond(bcondid,1:3) = { bcondid    bcondtype   bconddata}; 


    %================== GRAVITY LINE =======================
     elseif strfind(line,'*BFORCE')==1
         line = fgetl(fid);
         FE.bforce = str2num(line);
         
         
          %================== BOUNDARY LINE =======================
     elseif strfind(line,'*PFORCE')==1 % si encuentra la palabra despues del segundo caracter
         pforcename = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         pforcetype = regexp(line, '(?<=TYPE=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         pforceid=pforceid+1;
             line = fgetl(fid);
             idx = 0 ;

        while isempty(strfind(line,'*'))==1 
             idx = idx+1;
             pforcedata(idx,:) = str2num(line); %read section data and transform to number  
             curpos = ftell(fid); %get current position (to go back after the end of loop)
             line = fgetl(fid);
        end

        curpos = fseek(fid,curpos,'bof'); % move lo last file in order to avoid skeepint commands
        FE.pforce(pforceid,1:4) = { pforceid    pforcename   pforcetype   pforcedata};     
         
         
         
         
         

     %================== WIND LINE =======================
     elseif strfind(line,'*WIND')==1 % si encuentra la palabra despues del segundo caracter
         windid=windid+1; if windid>=2 ; disp('WARNING! Only one *WIND is allowed'); end
         windname = regexp(line, '(?<=NAME=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         bemflag  = regexp(line, '(?<=BEM=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
         windfname = regexp(line, '(?<=FILE=)\w*', 'match','once');

         if isempty(windfname)==0
             winddata = readiec (windfname); % Read .wnd file
         else
             line = fgetl(fid);
             winddata = str2num(line); %read section data and transform to number
         end

         FE.wind(1,:) = { windid    windname  bemflag   windfname}; 
         FE.wind(2,1:2) = { winddata };

         X.flag(10) = 1; 

             %================== PLOTS LINE =======================
     elseif strfind(line,'*PLOT')==1 % si encuentra la palabra despues del segundo caracter
        plotdata=[];
        plottype = regexp(line, '(?<=TYPE=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
        plotplane = regexp(line, '(?<=PLANE=)\w*', 'match','once');  % \w* matchea toda la palabra que acompaña
        curpos = ftell(fid); %get current position (to go back after the end of loop)
        line = fgetl(fid); 
        
        if isempty(strfind(line,'*'))==1
           while isempty(strfind(line,'*'))==1
             plotid=plotid+1;
             plotdata= str2num(line); % Read 4 parameters (OJO, solo 4) 
             curpos = ftell(fid); %get current position (to go back after the end of loop)
             FE.plots(plotid,1:3) = { plotid    plottype   plotdata}; 
             plotdata=[];
             line = fgetl(fid);           
           end
        else
           plotid=plotid+1;
           plotdata=0; % No variable
           FE.plots(plotid,1:4) = { plotid    plottype   plotdata  plotplane};
           curpos = ftell(fid); %get current position (to go back after the end of loop)
           line = fgetl(fid);
        end
           curpos = fseek(fid,curpos,'bof'); % move lo last file in order to avoid skeepint commands
      
     X.plotcount=plotid;
            
            
        %================== end of command list =======================
    end
    
    line = fgetl(fid);  %  READ NEW LINE
    
end 

X.secs=secid ; % Total Sections in model

%% ========================================================================
%                         PROPERTIES ASSIGNMENT
%==========================================================================

for i=1:X.secs
    
    if strcmp( FE.sec{2*i-1,3}, 'BEAM' ) == 1;       
        secmat=FE.sec(2*i-1,6);
        seclam=FE.sec(2*i-1,7);
        secmatid=strmatch(secmat,FE.mat(1:2:(2*matid),2),'exact');
        seclamid=strmatch(seclam,FE.lam(1:2:(2*lamid),2),'exact'); 
        FE.sec{2*i-1,6}=secmatid;
        FE.sec{2*i-1,7}=seclamid;
    end
    
    secset=FE.sec(2*i-1,4);
    secsetid=strmatch(secset,elset(:,1),'exact');
        
    if isempty(secsetid)==0
        [tf, eindex] = ismember(elset{secsetid,2}, E(:,1)); %#ok<ASGLU> % Hallo los indices de los elementos del set (por si estan desordenados
        E( eindex, 3 ) = i; % luego se altera si la seccion es variable
        E( eindex, 8 ) = i; % identificador del indice e la seccion
 
    else         
        FE.sec{2*i-1,4}={};  % Set to emtpy in order to avoid processing the section in pregen
        disp(strcat('Section ',' ',FE.sec(2*i-1,2), ' is not being used'))
    end
    
end    


fclose(fid);
