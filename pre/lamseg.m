% ======================================================================
%          CONSTITUTIVE MATRIX OF MULTIPLE MATERIALS LAMINATE
%
% 
% ======================================================================

function [CL,wseg,thseg] = lamseg(LS,MS,lamid,formulation)

CL = zeros (5,5); 
wseg = 0;
tsc = 0; % THICKNESS SHEAR COEFICIENT
%%=============================================================
%                  LAMINATE DATA
%%=============================================================
lyror  = LS{lamid,3}';      lyrs = length ( lyror);    % Orientatin and number of layers
lyrth = LS{lamid,4}';       thseg = sum(lyrth);         ninf = -thseg/2; % Thickness and Bottom of laminate
lyrmat = LS{lamid,5}';

% LAMINA INFO MATRIX (angle thickness zinf zsup)
LYR=zeros(lyrs,4);
LYR(:,1)= lyror;
LYR(:,2)= lyrth;
LYR(:,5)= lyrmat;
% Coordinates of every lamina
for lyr=1:lyrs
    nsup = ninf + LYR(lyr,2);
    LYR(lyr,3) = ninf ; %Z coordinate of bottom
    LYR(lyr,4) = nsup;
    ninf = nsup;
end


%%=============================================================
%             CONSTITUTIVE LOOP FOR EVERY LAMINA
%%=============================================================
for lyr = 1 : lyrs  

    %%=============================================================
    %                      MATERIAL DATA
    %%=============================================================
    mattype = MS{LYR(lyr,5),3};
    matprp  = MS{LYR(lyr,5),4};
    thlyr = LYR(lyr,4)-LYR(lyr,3);        % layer thickness
    % MATERIAL ISOTROPO
    if strcmp(mattype,'ISOTROPIC')==1 
        E1=matprp(1);   E2=E1;    
        v12=matprp(2); 
        G12=E1/(2+2*v12);   G13=G12;         G23=G12;
        wlyr = thlyr * matprp(3); % Weight per unit lenght of the layer
    end
    % LAMINA TRANSVERSALMENTE ISOTROPA CON EPT (G13=G12)
    if strcmp(mattype,'TRISOTROPIC')==1 
        E1=matprp(1);    E2=matprp(2);   
        G12=matprp(3);   G13=G12;        G23=matprp(4);
        v12=matprp(5); 
        wlyr = thlyr * matprp(6); % Weight per unit lenght of the layer
    end


    wseg = wlyr + wseg;  % Weight per unit length of laminate 
    
    mm=cosd( LYR(lyr,1) );     nn=sind( LYR(lyr,1) );  % directors   
     
    

    if strcmp(formulation,'BRB') == 1;

            % =====================================================================
            %       Constitutiva LOCAL de una lamina otrotropa  con EPT      
            % =====================================================================

            SL = [ 1/E1     -v12/E1    0        0       0       0   ;
                 -v12/E1     1/E2      0        0       0       0   ;
                   0           0       0        0       0       0   ;
                   0           0       0      1/G23     0       0   ;
                   0           0       0        0     1/G13     0   ;
                   0           0       0        0       0      1/G12 ];



            S11 = SL(1,1);      S22 = SL(2,2);     S12 = SL(1,2);
            S66 = SL(6,6);      S44 = SL(4,4);     S55 = SL(5,5);

            SG=zeros(6,6);

            SG(1,1) = S11*mm^4 + S22*nn^4 + (2*S12+S66)*nn^2*mm^2; % [Pa] = [N/m2]
            SG(1,2) = (S11+S22-S66)*nn^2*mm^2 + S12*(nn^4+mm^4);
            SG(2,2) = S11*nn^4 + S22*mm^4 + (2*S12+S66)*nn^2*mm^2;
            SG(1,6) = (2*S11-2*S12-S66)*nn*mm^3 - (2*S22-2*S12-S66)*mm*nn^3;
            SG(2,6) = (2*S11-2*S12-S66)*mm*nn^3 - (2*S22-2*S12-S66)*nn*mm^3;
            SG(6,6) = 2*(2*S11+2*S22-4*S12-S66)*mm^2*nn^2 + S66*(nn^4+mm^4);
            SG(4,4) = S55*nn^2 + S44*mm^2;
            SG(5,5) = S44*nn^2 + S55*mm^2;
            SG(4,5) = (S55-S44)*nn*mm;


        %     Ahora, eliminamos filas y columnas correspondientes a Ss Sns. ojo, hay EPT

%            SGEPT = [ SG(1,1)      0       SG(1,5) ;
%                         0       SG(4,4)      0    ;
%                      SG(1,5)      0       SG(5,5) ];
%                  
            SGEPT = [ SG(1,1)       0       SG(1,6) ;
                          0      SG(5,5)      0    ;
                      SG(1,6)      0        SG(6,6) ];
                  
                  

         % Invertimos para poner tensiones en funcion de deformaciones

            C =inv(SGEPT);


            A11 = C(1,1) * ( thlyr ) ;                                A16 = C(1,3) * ( thlyr );                                                     
            A61 = A16;                                                A66 = C(3,3) * ( thlyr ); 

            zk=(LYR(lyr,4)+LYR(lyr,3))/2;
            A55 = (5/4)*C(2,2) * ( thlyr - (4/thseg^2)*(thlyr*zk^2 + thlyr^3/12  )) ;             

            B11 = 0.5 * C(1,1) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );                B16 = 0.5 * C(1,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );
            B61 = B16;                                                           B66 = 0.5 * C(3,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );


            D11 = (1/3) * C(1,1) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );              D16 = (1/3) * C(1,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 ); 
            D61 = D16;                                                           D66 = (1/3) * C(3,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );


            CL =    [ A11   A16    0      B11    B16 ; 
                      A61   A66    0      B61    B66 ;
                       0     0    tsc*A55      0     0   ;
                      B11   B16    0      D11    D16 ;
                      B61   B66    0      D61    D66] + CL; % Laminate constitutive matrix



    end   % Barbero Formulation for Constants


    if strcmp(formulation,'SRV') == 1;       % SARAVIA Formulation for Constants    

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
        SG = inv(TE) * SL * TS; %#ok<MINV>


    %     Ahora, eliminamos filas y columnas correspondientes a Ss Sns. ojo, hay EPT

       SGEPT = [ SG(1,1)      0       SG(1,5) ;
                    0       SG(4,4)      0    ;
                 SG(1,5)      0       SG(5,5) ];

     % Invertimos para poner tensiones en funcion de deformaciones

      C =inv(SGEPT);



        A11 = C(1,1) * ( thlyr ) ;                                                        A16 = C(1,3) * ( thlyr );   
                                                       A55 = 0.25*C(2,2) * ( thlyr ) ;      
        A61 = A16;                                                                                        A66 = C(3,3) * ( thlyr ); 



        B11 = 0.5 * C(1,1) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );                B16 = 0.5 * C(1,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );
        B61 = B16;                                                           B66 = 0.5 * C(3,3) * ( LYR(lyr,4)^2 - LYR(lyr,3)^2 );


        D11 = (1/3) * C(1,1) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );              D16 = (1/3) * C(1,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 ); 
        D61 = D16;                                                           D66 = (1/3) * C(3,3) * ( LYR(lyr,4)^3 - LYR(lyr,3)^3 );



        CL =    [ A11   A16    0      B11    B16 ; 
                  A61   A66    0      B61    B66 ;
                   0     0   tsc*A55      0     0   ;
                  B11   B16    0      D11    D16 ;
                  B61   B66    0      D61    D66] + CL;

    end

end

