function IM=inercias(FE)

b=FE.sec{2};  h=FE.sec{3};  t=FE.sec{4};   



A = 2*(b + h)*t;   
Iy = (t*(h^3 + b*3*h^2 + b*t^2))/6; 
Iz = (t *(b^3 + 3 *b^2 * h + h *t^2))/6 ; 
J = ( (b + h)*t*((b + h)^2 + 4 *t^2))/6;

IM=[A Iy Iz J];