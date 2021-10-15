function xiT = xiT(V,phi)


% LINEALIZACION SEGUN MAKINEN, FUNCIONA BIEN DANDO MATRICES SIMETRICAS, ANTES LINEALIZABA SEGUN RITTO CORREA Y NO FUNCIONABA, NO ESTOY SEGURO SI EL PROBLEMA
% ESTA EN LA EXPRESION O BIEN EN QUE EXPRESAR LA LINEALIZACION EN FUNCION DE V =f_cross(se,a) Y NO EN FUNCION DE se COMO HABIA HECHO ANTES.
%   linealizacion de T 
%     V =f_cross(se,a);
% 
%     mphi = 1e-25 + sqrt(phi(1)^2+phi(2)^2+phi(3)^2);
%     a2=(1-cos(mphi))/(mphi^2);  
%     a3=(mphi-sin(mphi))/(mphi^3);
%     b1=(mphi*cos(mphi)-sin(mphi))/(mphi^3);     
%     b2=(mphi*sin(mphi)-2+2*cos(mphi))/(mphi^4);  
%     b3=(3*sin(mphi)-2*mphi-mphi*cos(mphi))/(mphi^5);
    % xiT =  a2*skew(V) + a3*(phi'*a)*dia3(1) + a3*(phi*a') + b1*(a*phi') - b2*((f_cross(phi,a))*phi') + b3*(phi'*a)*(phi*phi');
%     xiT = [ a3*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + b3*phi(1)^2*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + a3*a(1)*phi(1) + a(1)*b1*phi(1) + b2*phi(1)*(a(2)*phi(3) - a(3)*phi(2)),                                   a3*a(2)*phi(1) - V(3)*a2 + a(1)*b1*phi(2) + b2*phi(2)*(a(2)*phi(3) - a(3)*phi(2)) + b3*phi(1)*phi(2)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)),                                   V(2)*a2 + a3*a(3)*phi(1) + a(1)*b1*phi(3) + b2*phi(3)*(a(2)*phi(3) - a(3)*phi(2)) + b3*phi(1)*phi(3)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3));...
%                                    V(3)*a2 + a3*a(1)*phi(2) + a(2)*b1*phi(1) - b2*phi(1)*(a(1)*phi(3) - a(3)*phi(1)) + b3*phi(1)*phi(2)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)), a3*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + b3*phi(2)^2*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + a3*a(2)*phi(2) + a(2)*b1*phi(2) - b2*phi(2)*(a(1)*phi(3) - a(3)*phi(1)),                                   a3*a(3)*phi(2) - V(1)*a2 + a(2)*b1*phi(3) - b2*phi(3)*(a(1)*phi(3) - a(3)*phi(1)) + b3*phi(2)*phi(3)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3));...
%                                    a3*a(1)*phi(3) - V(2)*a2 + a(3)*b1*phi(1) + b2*phi(1)*(a(1)*phi(2) - a(2)*phi(1)) + b3*phi(1)*phi(3)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)),                                   V(1)*a2 + a3*a(2)*phi(3) + a(3)*b1*phi(2) + b2*phi(2)*(a(1)*phi(2) - a(2)*phi(1)) + b3*phi(2)*phi(3)*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)), a3*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + b3*phi(3)^2*(a(1)*phi(1) + a(2)*phi(2) + a(3)*phi(3)) + a3*a(3)*phi(3) + a(3)*b1*phi(3) + b2*phi(3)*(a(1)*phi(2) - a(2)*phi(1))];
%                                                 



mphi = 1e-25 + sqrt(phi(1)^2+phi(2)^2+phi(3)^2);
C1=(mphi*cos(mphi)-sin(mphi))/(mphi^3);
C2=(  mphi*sin(mphi)+2*cos(mphi)-2)/(mphi^4);
C3=(3*sin(mphi)-2*mphi-mphi*cos(mphi))/(mphi^5);
C4=(cos(mphi)-1)/(mphi^2);
C5=(mphi-sin(mphi))/(mphi^3); 

 Sphi = skew(phi);  SV=skew(V);

 xiT = C1*(V*phi') - C2*( Sphi*V )*phi' + C3 *( (phi'*V)*phi)*phi' - C4*SV + C5*( (phi'*V)*eye(3)+phi*V');