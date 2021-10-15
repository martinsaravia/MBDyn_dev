% ===================================================================
%    LAMINATE MATERIAL CONSTANTS FOR PLANE SHELL FORCE ASSUMTION   
%
% Assumtions: Equal thickness and even number of layers
% 
% Martin Saravia: 2011
% ===================================================================

function [C,rhon] = laminateEPS(FE,lamid,matid,thick)

%------------------------------------------------------------------------------------------------
matprp = FE.mat{matid*2,1};         lamprp = FE.lam{lamid*2,1};
mattype = FE.mat{matid*2-1,4};
%------------------------------------------------------------------------------------------------

e=thick; nlay=length(lamprp); ang=lamprp;


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


% _________________________________________________________________________
%                             RUTINA 1              
%              SALIDA: Ahij, Aij, Bij, Dij,  matrices con 
%                la SUMATORIA de cada Qij multipicado x
%            su cpte n/m/oplies (considera todo el laminado)
% _________________________________________________________________________

Aij = zeros(3,3);
Bij = zeros(3,3);
Dij = zeros(3,3);
Ahij = zeros(2,2);
for i = 1:nlay
    ang1 = ang(1,i);
    P = pi/180; 
%     eA = 1;
    mco = cos(P*ang1); nco=sin(P*ang1); 
    nu21 = v12*E2/E1;
% =========================================================================
% Coeficientes constitutivos para una LÁMINA ORTÓTROPA:: {Sg11}=[Qij]*{E1}
% =========================================================================
    Q11 = E1/(1-v12*nu21); % [Pa] = [N/m2]
    Q22 = E2/(1-v12*nu21);
    Q12 = nu21*E1/(1-v12*nu21);
    Q33 = G12;
    Q44 = G23;
    Q55 = G12;
% =======================================================================ok 
% Transformación a coordenadas curvilíneas:: {Sg} = [Qijr]*{E}
% =========================================================================
    Q11r = Q11*mco^4 + Q22*nco^4 + 2*(Q12+2*Q33)*nco^2*mco^2; % [Pa] = [N/m2]
    Q12r = (Q11+Q22-4*Q33)*nco^2*mco^2 + Q12*(nco^4+mco^4);
    Q22r = Q11*nco^4 + Q22*mco^4 + 2*(Q12+2*Q33)*nco^2*mco^2;
    Q13r = (Q11-Q12-2*Q33)*nco*mco^3 + (Q12-Q22+2*Q33)*mco*nco^3;
    Q23r = (Q11-Q12-2*Q33)*mco*nco^3 + (Q12-Q22+2*Q33)*nco*mco^3;
    Q33r = (Q11+Q22-2*Q12-2*Q33)*mco^2*nco^2 + Q33*(nco^4+mco^4);
    Q44a = Q55*nco^2 + Q44*mco^2;
    Q55a = Q44*nco^2 + Q55*mco^2;
    Q45a = (Q55-Q44)*nco*mco;
    Qrij = [Q11r,Q12r,Q13r;Q12r,Q22r,Q23r;Q13r,Q23r,Q33r];
    Qrij2 = [Q44a,Q45a;Q45a,Q55a];
% =======================================================================ok
    lplies = [-e/2:e/nlay:e/2];
    nplies = (lplies(i+1)-lplies(i));
    mplies = ((lplies(i+1)).^2-(lplies(i)).^2)/2;
    oplies = ((lplies(i+1)).^3-(lplies(i)).^3)/3;
     for k1 = 1:2
      for k2 = 1:2
      Ahij(k1,k2) = Ahij(k1,k2) + Qrij2(k1,k2)*nplies;
      end
     end
     for k3 = 1:3
      for k4 = 1:3
      Aij(k3,k4) = Aij(k3,k4) + Qrij(k3,k4)*nplies; % [Pa*m]
      Bij(k3,k4) = Bij(k3,k4) + Qrij(k3,k4)*mplies; % [Pa*m2]
      Dij(k3,k4) = Dij(k3,k4) + Qrij(k3,k4)*oplies; % [Pa*m3]
      end
     end       
end

CA= [Aij  Bij;Bij  Dij];

Delta = Aij(2,2)*Dij(2,2)-Bij(2,2)^2;
Del1 = (Bij(1,2)*Bij(2,2)-Aij(1,2)*Dij(2,2))/Delta;
Del2 = (Bij(2,2)*Bij(2,3)-Aij(2,3)*Dij(2,2))/Delta;
Del3 = (Dij(1,2)*Bij(2,2)-Bij(1,2)*Dij(2,2))/Delta;
Del4 = (Dij(2,3)*Bij(2,2)-Bij(2,3)*Dij(2,2))/Delta;
Del5 = -Ahij(1,2)/Ahij(1,1);
Del = [Del1,Del2,Del3,Del4,Del5];

% _________________________________________________________________________
%                              RUTINA 3
%                    Calculo de los A-ij, B-ij, D-ij; 
%            SALIDA: a11,a16,a33,b11,b16,b66,d11,d16,d66,rho,a55
% _________________________________________________________________________
K = -(Bij(2,2)^2)+Aij(2,2)*Dij(2,2);
Aijr = [Aij(1,1) + (-Aij(2,2)*Bij(1,2)^2 + 2*Aij(1,2)*Bij(1,2)*Bij(2,2) - Dij(2,2)*Aij(1,2)^2)/K, ...
        Aij(1,3) + (Aij(2,3)*Bij(1,2)*Bij(2,2) - Aij(2,2)*Bij(1,2)*Bij(2,3) + Aij(1,2)*Bij(2,2)*Bij(2,3) - Aij(1,2)*Aij(2,3)*Dij(2,2))/K ;
    Aij(1,3) + (Aij(2,3)*Bij(1,2)*Bij(2,2) - Aij(2,2)*Bij(1,2)*Bij(2,3) + Aij(1,2)*Bij(2,2)*Bij(2,3) - Aij(1,2)*Aij(2,3)*Dij(2,2))/K ,...
        Aij(3,3) + (2*Aij(2,3)*Bij(2,2)*Bij(2,3) - Aij(2,2)*Bij(2,3)^2 - Dij(2,2)*Aij(2,3)^2)/K];
Bijr = [Bij(1,1) + (Bij(2,2)*Bij(1,2)^2 - Aij(2,2)*Bij(1,2)*Dij(1,2) + Aij(1,2)*Bij(2,2)*Dij(1,2) - Aij(1,2)*Bij(1,2)*Dij(2,2))/K, ...
        Bij(1,3) + (Bij(1,2)*Bij(2,2)*Bij(2,3) - Aij(1,2)*Bij(2,3)*Dij(2,2) - Aij(2,2)*Bij(1,2)*Dij(2,3) + Aij(1,2)*Bij(2,2)*Dij(2,3))/K ;
    Bij(1,3) + (Bij(1,2)*Bij(2,2)*Bij(2,3) + Aij(2,3)*Bij(2,2)*Dij(1,2) - Aij(2,2)*Bij(2,3)*Dij(1,2) - Aij(2,3)*Bij(1,2)*Dij(2,2))/K, ...
        Bij(3,3) + (Bij(2,2)*Bij(2,3)^2 - Aij(2,3)*Bij(2,3)*Dij(2,2) + Aij(2,3)*Bij(2,2)*Dij(2,3) - Aij(2,2)*Bij(2,3)*Dij(2,3))/K];
Dijr = [Dij(1,1) + (2*Bij(1,2)*Bij(2,2)*Dij(1,2) - Aij(2,2)*Dij(1,2)^2 - Dij(2,2)*Bij(1,2)^2)/K, ...
        Dij(1,3) + (Bij(2,2)*Bij(2,3)*Dij(1,2) - Bij(1,2)*Bij(2,3)*Dij(2,2) + Bij(1,2)*Bij(2,2)*Dij(2,3) - Aij(2,2)*Dij(1,2)*Dij(2,3))/K ;
    Dij(1,3) + (Bij(2,2)*Bij(2,3)*Dij(1,2) - Bij(1,2)*Bij(2,3)*Dij(2,2) + Bij(1,2)*Bij(2,2)*Dij(2,3) - Aij(2,2)*Dij(1,2)*Dij(2,3))/K, ...
        Dij(3,3) + (2*Bij(2,2)*Bij(2,3)*Dij(2,3) - Dij(2,2)*Bij(2,3)^2 - Aij(2,2)*Dij(2,3)^2)/K];
A55r = Ahij(2,2) - Ahij(1,2)^2/Ahij(1,1);

a11 = Aijr(1,1);
a16 = Aijr(1,2);
a66 = Aijr(2,2);
a55 = A55r;  
b11 = Bijr(1,1); 
b16 = Bijr(1,2);
b66 = Bijr(2,2);
d11 = Dijr(1,1);
d16 = Dijr(1,2);
d66 = Dijr(2,2);

% A11 = Aij(1,1)-Aij(1,2)^2/Aij(2,2);
% A66 = Aij(3,3)-Aij(2,3)^2/Aij(2,2);
% D66 = Dij(3,3)-Dij(2,3)^2/Dij(2,2)

C=[  a11  a16   0   b11   b16 ; 
     a16  a66   0   b16   b66 ; 
     0     0   a55   0     0  ; 
     b11  b16   0   d11   d16 ; 
     b16  b66   0   d16   d66];


 