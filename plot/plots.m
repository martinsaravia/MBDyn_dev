%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%                             PLOTS
%                 
%                  MARTIN SARAVIA - 21/06/2010
%
%   Type 1: Plot 2D de una variable vs el tiempo
%   Type 2: Plot 3D de una deformada
%   Type 3: Animación.
%   
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


function plots(GUI,FE,X,E,N,U,Vel,Acc)


for plotid=1:X.plotcount
 
    plottype = FE.plots{plotid,2};
    pdataline = FE.plots{plotid,3};


    if pdataline(1,1)==1; var=U; name='Displacement'; end
    if pdataline(1,1)==2; var=Vel; name='Velocity'; end
    if pdataline(1,1)==3; var=Acc; name='Acceleration'; end
 
 
    if strcmp(plottype,'HISTORY')==1 % Historia Temporal
        dof=X.CTable(pdataline(1,2),pdataline(1,3));
        set(GUI.ploth(plotid),'xdata',X.T(2,:),'ydata',var(dof,1:X.tstep+1));
        drawnow
    end 


    if strcmp(plottype,'HF')==1
    %     forcetime=zeros(1,X.tstep);
    %     for i=1:200:X.tstep
    %       forcetime(1,i)=F(61,1)*X.T(3,i);
    %     end

        hold on

    %   vec=[11:40:4000]/dt;  vec=[1:10 vec]; vec=vec/40; % vec para Frame Invariant 

        vec=X.T(1,:);

        p1=plot( vec, var(dof,vec));
        set(p1, 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'Black');
        p2=plot( vec, var(dof+1,vec));
        set(p2, 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'Black');
        p3=plot( vec, var(dof+2,vec));
        set(p3, 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'Black');
    %     axis([ 0 100 -8 2.5]) %Frame Invariant
        drawnow

        set(gca,'FontSize',8,'FontName','Times','Box','off')
        xlabel('Revolutions','FontSize',8,'FontName','Times'); 
        ylabel('Displacements','FontSize',8,'FontName','Times'); 
        leg = legend('u','v','w','Location','NorthWest');
        set(leg,'Interpreter','none')

    end

    
    
    
    if strcmp(plottype,'SHAPE2')==1
      
        plotplane = FE.plots{plotid,4};

        ini=X.tstep; 

        for j=ini:X.tstep

            xx=N(X.actnodes,2)+U(X.CTable(X.actnodes,1),j);
            yy=N(X.actnodes,3)+U(X.CTable(X.actnodes,2),j);
            zz=N(X.actnodes,4)+U(X.CTable(X.actnodes,3),j);

           if strcmp(plotplane,'XY')==1; datax=xx; datay=yy; end
           if strcmp(plotplane,'XZ')==1; datax=xx; datay=zz; end
           if strcmp(plotplane,'YZ')==1; datax=yy; datay=zz; end

           maxcord = max([max(datax),max(datay)])*1.2; 
           scale=[-2 2];
           

           set(GUI.ploth(plotid),'xdata',datax,'ydata',datay,'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'Black');
           

        end
        set(GUI.axesh(plotid),'XLim',scale,'YLim',scale)
        xlabel(plotplane(1),'FontSize',8,'FontName','Times'); 
        ylabel(plotplane(2),'FontSize',8,'FontName','Times'); 


    end 

    if strcmp(plottype,'3D')==1
        subplot(2,2,sub,'Parent',GUI.uip2)
        set(gca,'FontSize',8,'FontName','Times','Box','off')
        axis equal

        if strcmp(type{2},'f')==1; ini=X.tstep; else ini=1; end

        for j=ini:X.tstep       
            actnodes=length(X.actnodes);
            xx=N(1:actnodes,2)+var(X.CTable(1:actnodes,1),j);
            yy=N(1:actnodes,3)+var(X.CTable(1:actnodes,2),j);
            zz=N(1:actnodes,4)+var(X.CTable(1:actnodes,3),j);
            p6=plot3(xx,yy,zz); set(p6, 'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'Black'); 
            axis([scale scale scale]); 
            title('Deformada 3D');  
    %         pause(0.00001)
        end
        xlabel('u','FontSize',8,'FontName','Times'); 
        ylabel('v','FontSize',8,'FontName','Times');
        zlabel('w','FontSize',8,'FontName','Times');
    end % Deformada 3D

    if strcmp (plottype,'SHAPE')==1

        % Figure Selection
      subplot(2,2,sub,'Parent',GUI.uip2)

    %       if sub<=4
    %         subplot(2,2,sub,'Parent',GUI.uip2)
    %       else
    %          figure(sub)
    %       end
    % 
    %       set(gca,'FontSize',8,'FontName','Times','Box','off')
    %       if strcmp(type{2},'f')==1; ini=X.tstep+1; fin=X.tstep+1; end  % Current Step Mode
    %       if strcmp(type{2},'a')==1; ini=1; fin=X.tstep+1; end          % Animation of full history
    %       if strcmp(type{2},'i')==1; ini=1; fin=1; end                  % Initial configuration
    %       if isnumeric(type{2})==1;  ini=type{2}; fin=type{2}; end                  % Initial configuration
    ini=X.tstep+1; 
    fin=X.tstep+1;

      for step=ini:fin

        idx=0;  
        for el=1:X.els       
            if E(el,2) == 3 || E(el,2) == 4 ||  E(el,2) == 5 || E(el,2) == 6   || E(el,2) == 8   || E(el,2) == 9
                idx = idx + 1;
                face(idx,:) = [E(el,4) E(el,5) E(el,5) E(el,4)];
            end
        end


        for nd =1:length(X.actnodes)
            nodo = X.actnodes(nd);
            udofs = X.CTable(nodo,1:3);
            vert(nodo,:) = N(nodo,2:4) + var(udofs,step)';
            mag(nodo,:)  = norm(var(udofs,step));
        end

    %  set (gca,'view',[60,20])  
       cla %clear current axes

    %        axis([scale scale scale])

        axis square
        view([-5,-5,5])

        patch('Vertices',vert,'Faces',face); %,'FaceVertexCData',mag,'FaceColor','interp','EdgeColor','interp'
        hold on
        xlabel('x','FontSize',8,'FontName','Times'); 
        ylabel('y','FontSize',8,'FontName','Times');
        zlabel('z','FontSize',8,'FontName','Times');    
    %     colorbar
    %     drawnow

      end

    end

    if strcmp (plottype,'MODO')==1

        modo = type{2}; amp=type{3};
        avec=zeros(X.sdof,1);  avec(:,1)=var(:,modo+1); %MODO+1 PORQUE LA PRIMER COLUMNA DE AV TIENE LAS FRECUENCIAS
        subplot(2,2,sub,'Parent',GUI.uip2);   set(gca,'FontSize',8,'FontName','Times','Box','off')       %     axis equal

    %     Nx=N(:,2);    Ny=N(:,3);    Nz=N(:,4);

        uu=avec(X.CTable(X.actnodes,1),1);
        vv=avec(X.CTable(X.actnodes,2),1);
        ww=avec(X.CTable(X.actnodes,3),1); 

        for t=1:0.2:10 
            amp2=amp*cos(t);
            xx=N(X.actnodes,2)+amp2*uu;
            yy=N(X.actnodes,3)+amp2*vv;
            zz=N(X.actnodes,4)+amp2*ww;

            URE=sqrt ( uu.^2 + vv.^2 + ww.^2 );   URE=[URE;URE];
            vert=[xx yy zz; xx yy zz+0.05];       face = [ E(:,4)   E(:,5)   E(:,5)+0.8    E(:,4)+0.8];

    %         set (gca,'view',[60,20])  
            cla %clear current axes
            axis([scale scale scale])
    %         axis([0 1 0 3 -1 1])
            view([5,5,5])
            patch('Vertices',vert,'Faces',face,'FaceVertexCData',URE,'FaceColor','flat'); %,'FaceVertexCData',URE,'FaceColor','interp'
            xlabel('x','FontSize',8,'FontName','Times'); 
            ylabel('y','FontSize',8,'FontName','Times');
            zlabel('z','FontSize',8,'FontName','Times');
            title(['Modo: ', num2str(modo),'  -  Freq(Hz): ', num2str(var(modo,1))],'FontSize',8,'FontName','Times')

    %         colorbar
            drawnow;  pause(0.1); 
        end

    end
end


% hold on
%  EXPORTAMOS LA FIGURA
% if strcmp(export{2},'YES')==1
%    filename=export{1};
%    export_fig ( (strcat(filename,'.eps')),-eps, figure(1)) 
% end



