function [T] = tangmap(phi)

phi_x = phi(1,1);       phi_y = phi(2,1);    phi_z = phi(3,1); 
penalty=1e-25;          mphi=penalty+sqrt(phi_x^2+phi_y^2+phi_z^2);
a1=sin(mphi)/mphi;      a2=(1-cos(mphi))/(mphi^2);      a3=(mphi-sin(mphi))/(mphi^3);

% Matriz T de Cardona (Ritto-Correa traspuesta)
T=[          a1 + a3*(phi_x)^2                   a3*(phi_x)*(phi_y) - a2*(phi_z)          a2*(phi_y) + a3*(phi_x)*(phi_z);  
       a3*(phi_x)*(phi_y) + a2*(phi_z)                 a1 + a3*(phi_y)^2                -(a2*(phi_x)) + a3*(phi_y)*(phi_z);
    -(a2*(phi_y)) + a3*(phi_x)*(phi_z)         a2*(phi_x) + a3*(phi_y)*(phi_z)                  a1 + a3*(phi_z)^2           ]';
















% tic
% % T de Makinen
%   TMak = a1* eye(3,3) - a2*skew(phi) + a3*(phi*phi')
%  toc
%  tic
%  % T de Cardona
%  TCard = eye(3,3) + ( (cos(mphi)-1)/mphi^2 )*skew(phi) + (1-(sin(mphi)/mphi))*( (skew(phi)*skew(phi))/mphi^2  ) 
%  toc
%  tic
%   % T de Ritto Correa 1
%  TRC1 = eye(3,3) + a2*skew(phi) + a3*skew(phi)*skew(phi)
%  toc
%  tic
%   % T de Rito Correa 2
%  TRC2 = a1* eye(3,3) + a2*skew(phi) + a3*(phi*phi')
%  toc
%  
% 
% identidad=skew(phi)*skew(phi)- ( (phi*phi')-(mphi^2)*eye(3,3) )






