function VelC = veldamp(X,FE,N,U,Vel)

VelC = zeros (length(Vel(:,1)),1);

omega= X.omega;


for nodo=1:53
    
    udofn=X.CTable(nodo,1:3);
    rdofn=X.CTable(nodo,4:6);
    
    x0 = N(nodo,2:4)' +  U(udofn,X.tstep+1);   
    
    vel = f_cross(omega,x0);
    
    VelC(udofn,1) = vel;
    VelC(rdofn,1) = omega;
    
end

VelC = Vel(:,X.tstep+1) - VelC;

% [I,C]=max(VelC)
% [I,C]=min(VelC)