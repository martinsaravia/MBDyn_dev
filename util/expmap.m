function [R] = expmap (phi)

% Given a finite rotation angle it computes the exponential map.

phi_x = phi(1);    phi_y = phi(2);    phi_z = phi(3); 
penalty=1e-25;       Mphi=penalty+sqrt(phi_x^2+phi_y^2+phi_z^2);



R = [1 + (2*(-(phi_y)^2 - (phi_z)^2)*sin(Mphi/2)^2)/Mphi^2 ...
(2*(phi_x)*(phi_y)*sin(Mphi/2)^2)/Mphi^2 - ((phi_z)*sin(Mphi))/Mphi ...
(2*(phi_x)*(phi_z)*sin(Mphi/2)^2)/Mphi^2 + ((phi_y)*sin(Mphi))/Mphi;

(2*(phi_x)*(phi_y)*sin(Mphi/2)^2)/Mphi^2 + ((phi_z)*sin(Mphi))/Mphi ...
1 + (2*(-(phi_x)^2 - (phi_z)^2)*sin(Mphi/2)^2)/Mphi^2 ...
(2*(phi_y)*(phi_z)*sin(Mphi/2)^2)/Mphi^2 - ((phi_x)*sin(Mphi))/Mphi;

(2*(phi_x)*(phi_z)*sin(Mphi/2)^2)/Mphi^2 - ((phi_y)*sin(Mphi))/Mphi ...
(2*(phi_y)*(phi_z)*sin(Mphi/2)^2)/Mphi^2 + ((phi_x)*sin(Mphi))/Mphi ...
1 + (2*(-(phi_x)^2 - (phi_y)^2)*sin(Mphi/2)^2)/Mphi^2];




