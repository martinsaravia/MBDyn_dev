disp ('---------- LAUNCHING PRE SCRIPT ----------')
FE.job =  'SNLWT100A2_NWP113';
FE.ext =  '.inps';
% FE.rundir = handles.dire;
atype = 'nonlinear'; %Choose the analysis type
handles=[]

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
%              plots(handles,FE,X,E,N,U,Vel,Acc); 
%               [X]=messages(handles,X,FE);
                disp([X.stpnfo(:,1) X.stpnfo(:,3)])

              break;    
             end
        end
       
     end
end

toc; disp (['Problem solved in: ' toc])
