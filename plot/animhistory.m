function animhistory(X,N,U)

%-------------------------------------------------
%          CONFIGURACION DE LA FIGURA
%-------------------------------------------------
h=figure(3);
set(h,'Position', [600, 200, 600, 400])
set(h, 'color', 'white')


%===========================================================
%              TIPO 1:  Animacion
%===========================================================

umax=max(max(U(1:6:X.sdof,:)));
vmax=max(max(U(2:6:X.sdof,:)));
wmax=max(max(U(3:6:X.sdof,:)));
umin=min(min(U(1:6:X.sdof,:)));
vmin=min(min(U(2:6:X.sdof,:)));
wmin=min(min(U(3:6:X.sdof,:)));

umax=max(umax); vmax=max(vmax); wmax=max(wmax);
umin=min(umin); vmin=min(vmin); wmin=min(wmin);
tol=1;

xmax=max(N(:,2))+tol; 
ymax=max(N(:,3))+tol; 
zmax=max(N(:,4))+tol;
xmin=min(N(:,2))-tol;
ymin=min(N(:,3))-tol;
zmin=min(N(:,4))-tol;

for i=1:length(U(1,:))
%       undef=plot(N(:,2),N(:,3)); %plot 2D undeformed shape
%       set(undef, 'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'Black');

 UX=N(:,2)+U(1:6:X.sdof,i);
 UY=N(:,3)+U(2:6:X.sdof,i);
 UZ=N(:,4)+U(3:6:X.sdof,i);


 plot3(UX,UY,UZ,'LineWidth',2)
%  plot(UX,UZ,'LineWidth',2)
 axis equal
%  axis([xmin+umin xmax+umax ymin+vmin ymax+vmax zmin+wmin zmax+wmax])
 drawnow
 hold on
 pause(0.1)   
end


