function [X]=buckle(X,N,K,KG,type,vals)


% ==========================================================================
%                             METODO DE  BATHE
%==========================================================================
if strcmp(type,'bathe')==1 
    
    K=0.5*(KG'+KG);
    KG=0.5*(KG'+KG);   % la simetrizo
    
    SYS=inv(-KG)*(K);
    options.disp=1; 
    [evec,eval]=eigs(SYS,vals,'SM',options);
    [X.evec2,X.eval2]=eig(K,KG);
    
    evec=real(evec); 
    eval=diag(real(eval));
    for i=1:vals
        X.eval(i,1)=(eval(i,1)-1)/(eval(i,1));
        X.evec(:,i)=bceliback(evec(:,i),X);
        X.eval2(i,1)=(X.eval2(i,1)-1)/(X.eval2(i,1));
    end    
     
end


% ==========================================================================
%                             METODO DE  CRISFIELD
%==========================================================================
if strcmp(type,'crisf')==1
       
    
%     KG=0.5*(KG'+KG);   % la simetrizo
%     K=0.5*(K'+K);
%     KM=K-KG;

    SYS=inv(-KG)*(K);

    options.disp=2; % options.tol=0.1;   %     options.p=5*vals;

    [evec,eval]=eigs(SYS,vals,'SM',options);
    evec=evec/(max(abs(evec)));
    evec=real(evec); X.eval=real(eval);
    eval
    %     [X.evec2,X.eval2]=eig(KM,-Kg); X.eval2=diag(X.eval2);


    for i=1:vals
    X.evec(X.adof,i)=evec(:,i);
    end   
end

% ==========================================================================
%                             METODO DE DIOS
%==========================================================================

if strcmp(type,'dios')==1
    [X.evec2,X.eval2]=eig(K,KG);

    X.eval2=diag(X.eval2);

    [X.eval2,k]=sort(diag(X.eval2));   % Ordeno los autovalores en orden ascendente.

    buckmode=X.evec2(:,k(8:vals+7));                               
    pcr=X.eval2; 

    % Loop para purgar numeros complejos y valores negativos
    temp=pcr;
    for i=1:length(temp)
      if real(temp(i))>=0.00001
      temp2(i)=real(temp(i));
      end
    end
    ind=find(temp2==0);
    temp2(ind)=[];
    X.pcr=temp2';

end


%     amp=1;
%     X.evec=amp*X.evec;
%     plot3(N(:,2),N(:,3),N(:,4)) %plot undeformed shape

%     umax=max(max(UT(1:6:X.sdof)));
%     vmax=max(max(UT(2:6:X.sdof)));
%     wmax=max(max(UT(3:6:X.sdof)));
%     umin=min(min(UT(1:6:X.sdof)));
%     vmin=min(min(UT(2:6:X.sdof)));
%     wmin=min(min(UT(3:6:X.sdof)));
%     umax=max(umax); vmax=max(vmax); wmax=max(wmax);
%     umin=min(umin); vmin=min(vmin); wmin=min(wmin);
% 
%     tol=1;
%     xmax=max(N(:,2))+tol; 
%     ymax=max(N(:,3))+tol; 
%     zmax=max(N(:,4))+tol;
%     xmin=min(N(:,2))-tol;
%     ymin=min(N(:,3))-tol;
%     zmin=min(N(:,4))-tol;



% for jjj=vals:vals;
%     UX=N(:,2)+X.evec(1:6:X.sdof,jjj);
%     UY=N(:,3)+X.evec(2:6:X.sdof,jjj);
%     UZ=N(:,4)+X.evec(3:6:X.sdof,jjj);
%     plot3(UX,UY,UZ,'LineWidth',1)
%     %    plot(UX,UZ)
% %         axis([xmin+umin xmax+umax ymin+vmin ymax+vmax zmin+wmin zmax+wmax])
%     hold on 
% end

