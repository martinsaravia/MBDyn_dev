function plotmodal(X,N,amp,modos,type)

if type==1
    %===========================================================
    %              TIPO 1:  Animacion MODO
    %===========================================================

    % plot3(N(:,2),N(:,3),N(:,4)) %plot undeformed shape
    umax=max(max(X.evec(1:6:X.sdof,:)));
    vmax=max(max(X.evec(2:6:X.sdof,:)));
    wmax=max(max(X.evec(3:6:X.sdof,:)));
    umin=min(min(X.evec(1:6:X.sdof,:)));
    vmin=min(min(X.evec(2:6:X.sdof,:)));
    wmin=min(min(X.evec(3:6:X.sdof,:)));
    umax=max(umax); vmax=max(vmax); wmax=max(wmax);
    umin=min(umin); vmin=min(vmin); wmin=min(wmin);
    tol=1;
    xmax=max(N(:,2))+tol; 
    ymax=max(N(:,3))+tol; 
    zmax=max(N(:,4))+tol;
    xmin=min(N(:,2))-tol;
    ymin=min(N(:,3))-tol;
    zmin=min(N(:,4))-tol;
    for t=1:0.2:10
        amp2=amp*cos(t);
        UX=N(:,2)+amp2*X.evec(1:6:X.sdof,modos);
        UY=N(:,3)+amp2*X.evec(2:6:X.sdof,modos);
        UZ=N(:,4)+amp2*X.evec(3:6:X.sdof,modos);
        plot3(UX,UY,UZ,'LineWidth',1)
        axis([xmin+umin xmax+umax ymin+vmin ymax+vmax zmin+wmin zmax+wmax])        
        drawnow;  pause(0.1); axis equal
    end
end


if type==2
    for modo=1:modos
        UX=N(:,2)+X.evec(1:6:X.sdof,modo);
        UY=N(:,3)+X.evec(2:6:X.sdof,modo);
        UZ=N(:,4)+X.evec(3:6:X.sdof,modo);
        plot3(UX,UY,UZ,'LineWidth',1)
        legend (num2str(X.eval))
        hold on; 
    end
end

