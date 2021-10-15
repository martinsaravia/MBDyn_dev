function spatch1d(X,E,N,var,sub,type,export,scale,name,view)
h=figure(1);
set(h,'Position', [200, 200, 600, 300])
set(h, 'color', 'white')
  
  set(gca,'FontSize',8,'FontName','Times','Box','off')
  if strcmp(type{2},'f')==1; ini=X.tstep+1; fin=X.tstep+1; end  % Current Step Mode
  if strcmp(type{2},'a')==1; ini=1; fin=X.tstep+1; end          % Animation of full history
  if strcmp(type{2},'i')==1; ini=1; fin=1; end                  % Initial configuration
  if isnumeric(type{2})==1;  ini=type{2}; fin=type{2}; end                  % Initial configuration
 
  for step=ini:fin

    idx=0;  
    for el=1:X.els       
        if E(el,2) == 3 || E(el,2) == 4        
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
    
%    axis(scale)
   axis equal

    patch('Vertices',vert,'Faces',face); %,'FaceVertexCData',mag,'FaceColor','interp','EdgeColor','interp'   
    xlabel('x','FontSize',8,'FontName','Times');  ylabel('y','FontSize',8,'FontName','Times');  zlabel('z','FontSize',8,'FontName','Times');    

%     colorbar
%     drawnow

  end