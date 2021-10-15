function [Mi] = esecvar (vec, e, phi)

% LINEALIZACION RITTO-CORREA 
Tinv=inv(tangmap(phi));
Ttinv=inv(tangmap(phi)');
Mi= Ttinv*(dtua(1,vec,e,phi))*Tinv + skew(vec)*skew(e); %#ok<MINV>

% LINEALIZACION DE GRUTTMAN ojo con el factor de penalty
phi_x = phi(1);    phi_y = phi(2);    phi_z = phi(3); 
penalty=1E-6;     mphi=penalty+sqrt(phi_x^2+phi_y^2+phi_z^2);

b=f_cross(e,vec);
c2=(mphi-sin(mphi))/(mphi^3);
c3=(mphi*sin(mphi)-2+2*cos(mphi))/(mphi^2*(cos(mphi)-1));
c6=(c3-c2)/mphi^2;
c9=-0.25+c3-0.25*(mphi^2)*c3^2;
c10=c2*(1-c9*mphi*mphi)*(b'*phi)-(e'*vec);
t=-c3*b + (c6+c2*c9)*(b'*phi)*phi;

Mi=0.5*(e*vec'+vec*e')+0.5*(t*phi'+phi*t')+c10*eye(3,3);


