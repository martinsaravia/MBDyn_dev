%% ==================================================
%  Producto Vectorial entre vectores a y b - a x b = cross
%
%          Todos los vectores son columna
%===================================================
function cross=f_cross(a,b)
% cross=zeros(3,1);
% cross(1) = a(2)*b(3) - a(3)*b(2);
% cross(2) = a(3)*b(1) - a(1)*b(3);
% cross(3) = a(1)*b(2) - a(2)*b(1);

cross=[a(2)*b(3)-a(3)*b(2) ; a(3)*b(1)-a(1)*b(3) ; a(1)*b(2)-a(2)*b(1)];