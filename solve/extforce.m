function F = extforce (FE,X,F,U,Ui,RT)



if isfield(FE, 'pforce')==1
    
    

DF = zeros(X.sdof,1);   MM = zeros(X.sdof,1);  MF = zeros(X.sdof,1);
    % % INCREMENTAL UPDATING
    % if X.iter==1 % Initial Iteration Assembly
    %     for i=1:length(FE.pforce(:,1))
    %         dofs = X.CTable(FE.pforce(i,1),1:6);     
    %         F(dofs,X.tstep+1) = F(dofs,X.tstep) + FE.pforce(i,2:7)';  %Fuerza TOTAL
    %     end
    % 
    % % ITERATIVE UPDATING
    % else 
    
    
% ----------------------------------------------------------------------
%        I N C R E M E N T A L     P O I N T     L O A D S 
% ---------------------------------------------------------------------- 
 for pforceid=1:length(FE.pforce(:,1))
     
    for i=1:length(FE.pforce{pforceid,4}(:,1))   % Loop por spatial applied moments.
        
        node = FE.pforce{pforceid,4}(i,1); 
        forcedata = FE.pforce{pforceid,4}(i,:); 
        
        nodeidx = node*3-2:node*3;  
        point=[0 0 0];  % Future support of excentric load
        
        udofs = X.CTable(node,1:3);      rdofs = X.CTable(node,4:6);
        
        force = forcedata(2:4)';       moment = forcedata(5:7)';
        
        % FORCES INCREMENTAL UPDATE
        DF(udofs,1) = force + DF(udofs,1);
        
        % MOMENTS ITERATIVE UPDATE
        TF = RT(nodeidx,1:3,2); % Incremental Tangmap at the node where force is applied    
        RF = RT(nodeidx,1:3,3); % Total Rotation at node
        
        % Applied Moments
        MM(rdofs,1) = TF * moment + MM(rdofs,1);
        
        
        % Moments of excentric load
%         if FE.pforce(i,11)==0;        % Fixed Force Flag
%             MF(rdofs,1) = TF * skew(point) * RF' * force + MF(rdofs,1);
%         elseif FE.pforce(i,11)==1;     %Follower Force Flag
%             MF(rdofs,1) = TF * skew(point) * force + MF(rdofs,1);
%         end

    end
 end
    
    F(:,X.tstep+1) = F(:,X.tstep) + DF + MM + MF;
    
else
    F(:,X.tstep+1) = F(:,X.tstep) ;
end

% Update of Total Force    



