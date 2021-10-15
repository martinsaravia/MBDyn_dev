% --------------------------------------------------
%         PLOT THE WARPING FUNCTION
% --------------------------------------------------

function plotd (GUI,NS,ES,wf)

% els = length(ES(:,1));
% nds = length(NS(:,1));

% sets = length (elset(:,1));


% for j=1:sets
%     setels = elset{j,2} ;  
%     for i=1:length (setels)  
%         el=setels(i);    n1=ES(el,2);   n2=ES(el,3);
%         edof = [DOF(n1,2) DOF(n2,2)];
%         p1=NS(n1,3:4);  p2=NS(n2,3:4);
%         line([0 0],[p1(1) p2(1)],[p1(2) p2(2)]) 
%         hold on
%         line([wf(edof(1),1) wf(edof(2),1)],[p1(1) p2(1)],[p1(2) p2(2)]) 
%     end
% end

% figure(2)
% NS
% wf
% 

% CROSS SECTION PLOT


psize = length(NS(:,3))+1;
datax=zeros(psize,1);
datay=zeros(psize,1);
dataw=zeros(psize,1);

datax (:,1)= [ NS(:,3) ; NS(1,3) ];
datay (:,1)= [ NS(:,4) ; NS(1,4) ];
dataw (:,1)= [ wf(:,1) ; wf(1,1) ];

axes(GUI.axes3); cla; 
line(datax,datay,'LineWidth',4) 
axis equal

axes(GUI.axes4); cla; 
line(datax,datay,dataw,'LineWidth',2) 
view(3)
axis equal



