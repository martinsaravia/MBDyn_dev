function plotblade(X,FE,E,N,id)

didx=id*2;

foilshape = FE.sec{didx,2};     bladeshape = FE.sec{didx,3};

long = max(FE.sec{didx,3}(:,1));

laxis=0:1:long;

nseg = length(laxis);

chord = interp1( bladeshape(:,1) , bladeshape(:,3) , laxis, 'cubic' );

points = length(foilshape(:,1));

for i=1:nseg
    
    xsec = chord(i) * foilshape;
    
    zeta = laxis(i) * ones(points,1);
    
    
    plot3 ( xsec(:,1), xsec(:,2), zeta)
    axis equal
    hold on
    
end