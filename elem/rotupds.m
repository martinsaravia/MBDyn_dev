function [RT] = rotupds(X,RT,Ui,fase)


if strcmp(fase,'corr')==1

    for i = 1:length(X.actnodes)
        node = X.actnodes(i);
        idx = node*3-2 : node*3;  
        rdof = X.CTable(node,4:6);

        Dphi = Ui(rdof,X.tstep+1); 

        RT(idx,1:3,1) = expmap(Dphi);  % Incremental Rotation Matrix
        RT(idx,1:3,2) = tangmap(Dphi); % Incremental Tangential Transforation
        RT(idx,1:3,3) = RT(idx,1:3,4) * RT(idx,1:3,1); % Total Rotation = Rref*Rinc
        
    end

    % Reference Rotation Update
    if X.conv==1
        for i = 1:length(X.actnodes)
            node = X.actnodes(i);
            idx = node*3-2 : node*3;         
            RT(idx,1:3,1) = eye(3);
            RT(idx,1:3,2) = eye(3);  % No la inicializo a cero porqeu la usa en advance para calcular las velocidades angulares iniciales
            RT(idx,1:3,4) = RT(idx,1:3,3); % Update Reference Rotation Matrix
        end
    end
    
end

if strcmp(fase,'pred')==1  % Rotation update for the prediction
    for i = 1:length(X.actnodes)
        node = X.actnodes(i);
        idx = node*3-2 : node*3;  
        rdof = X.CTable(node,4:6);

        Dphi = Ui(rdof,X.tstep+1); 

        RT(idx,1:3,1) = expmap(Dphi);  % Incremental Rotation Matrix
        RT(idx,1:3,2) = tangmap(Dphi); % Incremental Tangential Transforation
        RT(idx,1:3,3) = RT(idx,1:3,4) * RT(idx,1:3,1); % Total Rotation = Rref*Rinc
        
    end
end