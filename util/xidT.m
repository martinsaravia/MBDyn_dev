function xidT = xidT(dphi,phi)


% LINEALIZACIÓN SEGÚN MAKINEN                                        

mphi = 1e-25 + sqrt(phi(1)^2+phi(2)^2+phi(3)^2);
C1=(mphi*cos(mphi)-sin(mphi))/(mphi^3);
C2=(  mphi*sin(mphi)+2*cos(mphi)-2)/(mphi^4);
C3=(3*sin(mphi)-2*mphi-mphi*cos(mphi))/(mphi^5);
C4=(cos(mphi)-1)/(mphi^2);
C5=(mphi-sin(mphi))/(mphi^3); 

Sphi = skew(phi);  Sdphi=skew(dphi);  I = eye(3);

xidT = C1*(phi'*dphi)*I - C2*( phi'*dphi )*Sphi + C3 *((phi'*dphi)*phi)*phi' + C4*Sdphi + C5*(dphi*phi' + phi*dphi');