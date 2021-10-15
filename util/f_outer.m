%%==============================================================================================
%
%                          CALCULO DE TENSOR DE PRODUCTO EXTERNO
%    IN:  vectores
%    OUT: tensor de producto externo
%==============================================================================================


function outer=f_outer(v1,v2)

 
outer = [ v1(1)*v2(1)     v1(1)*v2(2)    v1(1)*v2(3);
          v1(2)*v2(1)     v1(2)*v2(2)    v1(2)*v2(3);
          v1(3)*v2(1)     v1(3)*v2(2)    v1(3)*v2(3)];    
      
end