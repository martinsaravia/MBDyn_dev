function [sw] = int1d (nips)

if nips == 1 
    
    sw(1,1) = 0.0;
    sw(2,1) = 2.0;

elseif nips == 2 

    sw(1,1) = -sqrt(1/3);                  
    sw(1,2) =  sqrt(1/3);
    sw(2,1) = 1.0;
    sw(2,2) = 1.0;

else
    disp('Too much integration points - extend routine')
end