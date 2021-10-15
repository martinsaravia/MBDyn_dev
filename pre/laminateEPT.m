% ===================================================================
%    LAMINATE MATERIAL CONSTANTS FOR PLANE STRESS ASSUMTION   
%
% Assumtions: Equal thickness and even number of layers
% 
% Martin Saravia: 2011
% ===================================================================

function [CL,rhon] = laminateEPT(FE, lamid, matid, thick)


%------------------------------------------------------------------------------------------------
matprp = FE.mat{matid*2,1};         lamprp = FE.lam{lamid*2,1};
mattype = FE.mat{matid*2-1,4};
%------------------------------------------------------------------------------------------------


CL = zeros (5,5); 
lyrs=length(lamprp);             % Number of layers
lyrth  = thick / lyrs;  % Reparto espesores iguales para cada lamina
ninf = -thick/2; %Bottom of laminate


% LAMINA INFO MATRIX (angle thickness zinf zsup)
LYR=zeros(lyrs,4);
LYR(:,1)= lamprp';
LYR(:,2)= lyrth;

for lyr=1:lyrs
    nsup = ninf + lyrth;
    LYR(lyr,3) = ninf ; %Z coordinate of bottom
    LYR(lyr,4) = nsup;
    ninf = nsup;
end

%%=============================================================
%                 CONSTITUTIVE LAW SELECTION
%%=============================================================
% MATERIAL ISOTROPO
if strcmp(mattype,'ISOTROPIC')==1 
    E1=matprp(1);   E2=E1;    
    v12=matprp(2); 
    G12=E1/(2+2*v12);   G13=G12;         G23=G12;
    rhon = thick * matprp(3); %Density per unit thickness
end
% LAMINA TRANSVERSALMENTE ISOTROPA CON EPT (G13=G12)
if strcmp(mattype,'TRISOTROPIC')==1 
    E1=matprp(1);    E2=matprp(2);   
    G12=matprp(3);   G13=G12;        G23=matprp(4);
    v12=matprp(5); 
    rhon = thick * matprp(6); %Density per unit thickness
end

%%=============================================================
%                  LAMINATE STIFFNESS
%%=============================================================

for lyr = 1 : lyrs   %LOOP FOR EVERY LAMINA
    
    mm=cosd( LYR(lyr,1) );     nn=sind( LYR(lyr,1) );  % directors      
                       
    % =====================================================================
    %       Constitutiva LOCAL de una lamina otrotropa  con EPT      
    % =====================================================================
    
    SL = [ 1/E1     -v12/E1        0       0       0   ;
         -v12/E1     1/E2          0       0       0   ;
           0           0         1/G23     0       0   ;
           0           0           0     1/G13     0   ;
           0           0           0      0      1/G12 ];
       

    % =====================================================================
    %                     TRANSFORMACIONES      
    % =====================================================================  
    % Transforamcion de tensiones
    TS = [ mm^2   nn^2          0     0     2*mm*nn    ;
           nn^2   mm^2          0     0    -2*mm*nn    ;
           0       0            mm   -nn       0       ;
           0       0            nn    mm       0       ;
        -mm*nn   mm*nn          0     0     mm^2-nn^2   ];
    
    TE = [ mm^2   nn^2          0     0     mm*nn    ;
           nn^2   mm^2          0     0    -mm*nn    ;
           0       0            mm   -nn     0       ;
           0       0            nn    mm     0       ;
        -2*mm*nn 2*mm*nn        0     0    mm^2-nn^2 ];


    % TRANSFORMO LA CONSTITUTIVA AL SISTEMA GLOBAL
    SG = inv(TE) * SL * TS;
    
    
    
%     Ahora, eliminamos filas y columnas correspondientes a Ss Sns. ojo, hay EPT
    
   SGEPT = [ SG(1,1)      0       SG(1,5) ;
                0       SG(4,4)      0    ;
             SG(1,5)      0       SG(5,5) ];
      
 % Invertimos para poner tensiones en funcion de deformaciones
 
  C =inv(SGEPT);

    
    A11 = C(1,1) * ( LYR(lyr,4)-LYR(lyr,3) ) ;                                                        A16 = C(1,3) * ( LYR(lyr,4)-LYR(lyr,3) );   
                                                   A55 = C(2,2) * ( LYR(lyr,4)-LYR(lyr,3) ) ;      
    A61 = A16;                                                                                        A66 = C(3,3) * ( LYR(lyr,4)-LYR(lyr,3) ); 
                
                
    
    B11 = 0.5 * C(1,1) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );                B16 = 0.5 * C(1,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );
    B61 = B16;                                                           B66 = 0.5 * C(3,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );
    
    
    D11 = (1/3) * C(1,1) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );              D16 = (1/3) * C(1,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 ); 
    D61 = D16;                                                           D66 = (1/3) * C(3,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );
                
                

    CL =    [ A11   A16    0      B11    B16 ; 
              A61   A66    0      B61    B66 ;
               0     0    A55      0     0   ;
              B11   B16    0      D11    D16 ;
              B61   B66    0      D61    D66] + CL;
           

end



