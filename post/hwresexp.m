%% 
%       EXPORT HW ASCII RESULTS FILE

function  hwresexp (X,FE,E,N,U,Vel,Acc,AV,type)



resfile = fopen( [FE.job '.hwascii'],'w');


fprintf(resfile,'ALTAIR ASCII FILE\n');
fprintf(resfile,'$DELIMITER =    <SPACE>/<TAB>/\n');
fprintf(resfile,'$TITLE = %s\n',FE.job);

if strcmp(type,'dyn')==1
    %==================== EXPORT DISPLACEMENTS ========================;
    fprintf(resfile,'$BINDING = NODE\n');
    fprintf(resfile,'$ID_POOL = Beam\n');
    fprintf(resfile,'$COLUMN_INFO = ENTITY_ID \n');
    fprintf(resfile,'$RESULT_TYPE =  Displacement(v), Rotation(v)  \n');
    for t=1:X.tstep+1    
        time = X.T(2,t);  nodes = length(X.actnodes);  
        fprintf(resfile,'$TIME =  %.6f \n', time);
        table = zeros (nodes,7);
        for inode=1:X.ands
            node=X.actnodes(inode);
            dofs = X.CTable(node,1:6); 
            table(inode,1)=node;
            table(inode,2:7) = U(dofs,t); 
        end 
        table = table';
        fprintf(resfile,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',table);
    end

    if strcmp (FE.step{1,3},'DYNAMIC') ==1
    %==================== EXPORT VELOCITIES ========================;
        fprintf(resfile,'$BINDING = NODE\n');
        fprintf(resfile,'$ID_POOL = Beam\n');
        fprintf(resfile,'$COLUMN_INFO = ENTITY_ID \n');
        fprintf(resfile,'$RESULT_TYPE =  VelocityT(v), VelocityR(v)  \n');
        for t=1:X.tstep+1    
            time = X.T(2,t);  nodes = length(X.actnodes);  
            fprintf(resfile,'$TIME =  %.6f \n', time);
            table = zeros (nodes,7);
            for inode=1:X.ands
                node=X.actnodes(inode);
                dofs = X.CTable(node,1:6); 
                table(inode,1)=node;
                table(inode,2:7) = Vel(dofs,t); 
            end 
            table = table';
            fprintf(resfile,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',table);
        end
        % 
        %==================== EXPORT ACCELERATIONS ========================;
        fprintf(resfile,'$BINDING = NODE\n');
        fprintf(resfile,'$ID_POOL = Beam\n');
        fprintf(resfile,'$COLUMN_INFO = ENTITY_ID \n');
        fprintf(resfile,'$RESULT_TYPE =  AccelerationT(v), AccelerationR(v)  \n');
        for t=1:X.tstep+1    
            time = X.T(2,t);  nodes = length(X.actnodes);  
            fprintf(resfile,'$TIME =  %.6f \n', time);
            table = zeros (nodes,7);
            for inode=1:X.ands
                node=X.actnodes(inode);
                dofs = X.CTable(node,1:6); 
                table(inode,1)=node;
                table(inode,2:7) = Acc(dofs,t); 
            end 
            table = table';
            fprintf(resfile,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',table);
        end
    end
    %==================== EXPORT TRIADS ========================;
    fprintf(resfile,'$BINDING = ELEMENT\n');
    fprintf(resfile,'$ID_POOL = Beam\n');
    fprintf(resfile,'$COLUMN_INFO = ENTITY_ID \n');
    fprintf(resfile,'$RESULT_TYPE =  Triad_x(v), Triad_y(v), Triad_z(v)  \n');
    
    for t=1:X.tstep+1    
        time = X.T(2,t);  
        fprintf(resfile,'$TIME =  %.6f \n', time);
        table = zeros (X.els,10);
        for el=1:X.els
            table(el,1) = E(el,1);
            idxs=el*3-2:el*3; 
            table(el,2:4) = AV(idxs,1,t)';
            table(el,5:7) = AV(idxs,2,t)';
            table(el,8:10) = AV(idxs,3,t)'; 
        end
        table = table';
        fprintf(resfile,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',table);
    end


end



% fprintf(resfile,'$BINDING = ELEMENT\n');
% fprintf(resfile,'$ID_POOL = Beam\n');
% fprintf(resfile,'$COLUMN_INFO = ENTITY_ID \n');
% 
% fprintf(resfile,'# ------- TRIADS RESULTS -------\n');
% for t=1:X.tstep    
%     time = X.T(2,t);  nodes = length(X.actnodes);  dofs = X.sfdof ;
%     fprintf(resfile,'$RESULT_TYPE =  Triad_x(v), Triad_y(v), Triad_z(v)  \n');
%     fprintf(resfile,'$TIME =  %.6f \n', time);
%     table = zeros (nodes,10);
%     table(:,1) = X.actnodes; % Node IDs
%     length(table(:,2) )
%     length( U(1:6:dofs,t) )
%     table(:,2) = U(1:6:dofs,t);    table(:,3) = U(2:6:dofs,t);    table(:,4) = U(3:6:dofs,t); 
%     table(:,5) = U(4:6:dofs,t);    table(:,6) = U(5:6:dofs,t);    table(:,7) = U(6:6:dofs,t); 
%     table = table';
%     fprintf(resfile,'%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n',table);
% end
% fclose(resfile);

% fprintf(resfile,' \n');
% fprintf(resfile,' \n');
% fprintf(resfile,' \n');
% fprintf(resfile,' \n');
% fprintf(resfile,' \n');
% fprintf(resfile,' \n');
% 
% $DELIMITER=    <SPACE>/<TAB>/,    (except $ and #)
% $TITLE=
% $SUBCASE_ID =           id        “label”
% $BINDING = ELEMENT or NODE
% $LAYER_INFO =        num layers        “layer label”
% $ID_POOL= “Nodes” or “Elements”
% $SYS_ID = -1/0/user defined
% $COLUMN_INFO=        ENTITY_ID GRID_ID
% $RESULT_TYPE =  RES(v)(R/I)(Anim Source), RES(t)(M/P), RES(s)        
% $FREQUENCY =  1.00 Hz
% $MODE = 1
% $TIME= 10.0 sec
% Data blocks OR
% $INCLUDE  file path (which has the data)