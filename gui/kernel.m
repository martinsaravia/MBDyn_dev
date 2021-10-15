
clc;   clear -handles
disp ('---------- LAUNCHING PRE SCRIPT ----------')
FE.job =  get(handles.fileedit,'String');
FE.ext =  handles.exte;
FE.rundir = handles.dire;


%%----------------------------------------------------------------------------------
%                      CHOOSE ANALYSIS TYPE
%----------------------------------------------------------------------------------

if strcmp(handles.exte ,'.sinp')==1
    atype = 'sectional'; %Choose the analysis type
elseif strcmp(handles.exte ,'.inps')==1
    atype = 'nonlinear'; %Choose the analysis type
end

%%----------------------------------------------------------------------------------
%                      SECTIONAL ANALYSIS ONLY
%----------------------------------------------------------------------------------
if strcmp(atype,'sectional')==1
 filename=[FE.job FE.ext];   
 tic; disp (['Processing: ' filename])
 [NS,ES,D,DM,LS,MS,wf,info] = sectgen(filename, 1, 1,'outer','NO', 'LINEAR');
 plotd (handles,NS,ES,wf);
 toc; disp (['Sectional Properties solved in: ' toc])
end



%%----------------------------------------------------------------------------------
%                    DYNAMIC ANALYSIS
%----------------------------------------------------------------------------------
if strcmp(atype,'dynamic')==1
    
tic; disp (['Processing: ' FE.job FE.ext])
[FE,X,N,E,D,DM,WF,F,FD,FI,FM,FA,SE,EE,U,Ui,DU,UT,Vel,Acc,Ome,Gam,A0,AV,RT]=pregen(FE);
toc; disp (['Pre done in: ' toc])
     
X.tn =  AV;


tic; disp ('---------- LAUNCHING KERNEL ----------')
if FE.step{2,2}(3) == 1  % NONLINEAR GEOMETRIC ANALYSIS
    disp ('Nonlinear analysis has started...')
      for tstep = 1:X.totsteps; 
        X.tstep=tstep;
        [F,X,U,Ui,Vel,Acc,Ome,Gam,A0,RT]=advance(X,FE,F,U,Ui,Vel,Acc,Ome,Gam,A0,RT); 
        for iter = 1:FE.param(4); X.iter=iter;
            [X,M,C,K,KG,F,FD,FI,FM,FA,SE,EE,AV]=tangent(FE,X,E,N,U,Ui,Vel,Acc,Ome,Gam,DU,D,DM,AV,F,FD,FI,FM,FA,SE,EE,RT);
            [X,U,Ui,Vel,Acc,Ome,Gam,RT] = iterate (FE,X,N,U,Ui,Vel,Acc,Ome,Gam,K,C,M,F,FD,FI,FM,FA,RT);
             if X.conv==1; 
                 plots(handles,FE,X,E,N,U,Vel,Acc); 
                 [X]=messages(handles,X,FE);
              break;    
             end
        end
       if get(handles.stopbutton,'UserData') == 1;  break;   end
       
     end
end

toc; disp (['Problem solved in: ' toc])
end



%%----------------------------------------------------------------------------------
%                   STATIC ANALYSIS
%----------------------------------------------------------------------------------
if strcmp(atype,'nonlinear')==1
    
    disp (['Processing: ' FE.job FE.ext])
   [FE,X,N,E,D,DM,WF,F,FD,FI,FM,FA,SE,EE,U,Ui,DU,UT,Vel,Acc,Ome,Gam,A0,AV,RT]=pregen(FE);
   disp ('Pre Succesfully Completed...')

    disp ('---------- LAUNCHING KERNEL ----------')
    if FE.step{2,2}(3) == 1  % NONLINEAR GEOMETRIC ANALYSIS
        disp ('Nonlinear analysis has started...')
        tic
          for tstep = 1:X.totsteps; 
            X.tstep=tstep;
            [F,X,U,Ui,Vel,Acc,Ome,Gam,A0,RT]=advance(X,FE,F,U,Ui,Vel,Acc,Ome,Gam,A0,RT); 
            for iter = 1:FE.param(4); X.iter=iter;
                 [X,M,C,K,KG,F,FD,FI,FM,FA,SE,EE,AV]=tangent(FE,X,E,N,U,Ui,Vel,Acc,Ome,Gam,DU,D,DM,AV,F,FD,FI,FM,FA,SE,EE,RT);
                 [X,U,Ui,Vel,Acc,Ome,Gam,RT] = iterate (FE,X,N,U,Ui,Vel,Acc,Ome,Gam,K,C,M,F,FD,FI,FM,FA,RT);
                 if X.conv==1;  plots(handles,FE,X,E,N,U,Vel,Acc);  [X]=messages(handles,X,FE); break;    end               
            end
            if get(handles.stopbutton,'UserData') == 1;  break;   end
          end
       toc
    end
disp ('Analysis Succesfully Completed...')
end


save(FE.job)



% if strcmp(atype,'linear')==1
%     disp (' Linear Analysis has Started')
%     [X,ML,KL,KG,F,FD,FI,FM,FA,SE,EE,AV]=tangent(FE,X,E,N,U,Ui,Vel,Acc,DU,D,DM,AV,F,FD,FI,FM,FA,SE,EE,RT);  
%     for tstep = 1:tsteps; X.tstep=tstep;
%         for iter = 1:1; X.iter=iter;       
%         [F,X,U,Ui,Vel,Acc,A0,RT]=advance(X,FE,F,U,Ui,Vel,Acc,A0,RT); 
%         [X,M,K,KG,F,FD,FI,FM,FA,SE,EE,AV]=tangent(FE,X,E,N,U,Ui,Vel,Acc,DU,D,DM,AV,F,FD,FI,FM,FA,SE,EE,RT);     
%         [X,U,Ui,Vel,Acc,Ome,Gam,RT] = iterate (FE,X,N,U,Ui,Vel,Acc,K,C,M,F,FD,FI,FM,FA,RT);  
%              if X.conv==1; 
%                           plots(X,E,N,U,  1,{'HI' 64 },{'' ''},gscale,'Disp'); 
%                           plots(X,E,N,U,  2,{'2D' 'f' 'u' 'w' },{'' ''},gscale,'Disp');
%                           plots(X,E,N,Acc,3,{'HI' 64 },{'' ''},gscale,'Acc');
%                           plots(X,E,N,U,  4,{'PATCH' 'f'},{'' ''},gscale);... 
%                           
%              [X] = messages(X,FE); 
%              break;    
%              end
%         end
%     end
% 
% end
% 
% 
% disp ('Analysis Succesfully Completed...')


