function plot_rollup(X,N,UT,F,amp,type,dof)

%-------------------------------------------------
%          CONFIGURACION DE LA FIGURA
%-------------------------------------------------
h=figure(1);
set(h,'Position', [600, 200, 600, 400])
set(h, 'color', 'white')


%===========================================================
%              TIPO 1:  Plot 1D
%===========================================================
if type==1
  U=UT;
  forcetime=zeros(1,X.tstep);
  for i=1:X.tstep
      forcetime(1,i)=F(61,1)*X.T(3,i);
  end
  hold on
  plot(-U(dof,1:X.tstep),X.T(3,1:X.tstep))
  drawnow
end


%===========================================================
%              TIPO 1:  Animacion
%===========================================================
if type==2
  
  UT=amp*UT;

  totalsteps=size(UT); 

  if length(totalsteps)==3
    totalsteps=totalsteps(3); %el array 3D aparece solo para mas de un load case
  else
   totalsteps=1;
  end

  umax=zeros(totalsteps,1); vmax=zeros(totalsteps,1); wmax=zeros(totalsteps,1); 
  umin=zeros(totalsteps,1); vmin=zeros(totalsteps,1); wmin=zeros(totalsteps,1);

  if totalsteps>=2
    for i=1:totalsteps
      umax(i)=max(max(UT(1:6:X.sdof,:,i)));
      vmax(i)=max(max(UT(2:6:X.sdof,:,i)));
      wmax(i)=max(max(UT(3:6:X.sdof,:,i)));
      umin(i)=min(min(UT(1:6:X.sdof,:,i)));
      vmin(i)=min(min(UT(2:6:X.sdof,:,i)));
      wmin(i)=min(min(UT(3:6:X.sdof,:,i)));
    end
  end
  if totalsteps==1
      umax=max(max(UT(1:6:X.sdof,:)));
      vmax=max(max(UT(2:6:X.sdof,:)));
      wmax=max(max(UT(3:6:X.sdof,:)));
      umin=min(min(UT(1:6:X.sdof,:)));
      vmin=min(min(UT(2:6:X.sdof,:)));
      wmin=min(min(UT(3:6:X.sdof,:)));
  end

  umax=max(umax); vmax=max(vmax); wmax=max(wmax);
  umin=min(umin); vmin=min(vmin); wmin=min(wmin);
  tol=1;
  xmax=max(N(:,2))+tol; 
  ymax=max(N(:,3))+tol; 
  zmax=max(N(:,4))+tol;
  xmin=min(N(:,2))-tol;
  ymin=min(N(:,3))-tol;
  zmin=min(N(:,4))-tol;


    
%       undef=plot(N(:,2),N(:,3)); %plot 2D undeformed shape
%       set(undef, 'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'Black');
     hold  on
     

                  
%    subplot(1,2,1); 
%    plot3(UX,UY,UZ,'LineWidth',2)
%    axis([xmin+umin xmax+umax ymin+vmin ymax+vmax zmin+wmin zmax+wmax])
%    drawnow
%    pause(0.05)   
%    

% %    subplot(1,2,1); 
%      tstep=5;
%      UX=N(:,2)+UT(1:6:X.sdof,tstep);
%      UY=N(:,3)+UT(2:6:X.sdof,tstep);
%      UZ=N(:,4)+UT(3:6:X.sdof,tstep);
%      def1=plot(UX,-UZ);  axis equal
%      set(def1, 'LineStyle', '--', 'LineWidth', 2.0, 'Color', 'Black');
% %      set(def1, 'Marker', 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'MarkerSize',1.0);
%      drawnow
% %      pause(0.05)

hold on

     tstep=X.tstep;
     UX=N(:,2)+UT(1:6:X.sdof,tstep);
     UY=N(:,3)+UT(2:6:X.sdof,tstep);
     UZ=N(:,4)+UT(3:6:X.sdof,tstep);
     def1=plot(UX,-UZ);  axis equal
     set(def1, 'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'Black');
%      set(def1, 'Marker', 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'MarkerSize',1.0);
     drawnow
%      pause(0.05)
  
    axis([-25 25 0 35])
%-------------------------------------------------
%          CONFIGURACION DE LA FIGURA
%-------------------------------------------------
    set(gca,'FontSize',8,'FontName','Times','Box','on')
    xlabel('x\it','FontSize',8,'FontName','Times'); 
    ylabel('z\it','FontSize',8,'FontName','Times'); 
    leg = legend('M1','M2','Location','NorthWest');
%     set(leg,'Interpreter','none')
   
  end


%  filename='rollup';
%  export_fig ( (strcat(filename,'.eps')),-eps, h) 
