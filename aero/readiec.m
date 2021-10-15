%% ====================================================
%               READ wnd files from IECWind
% ====================================================
function VWT = readiec (fname)

VWT = zeros(1,8);
fid = fopen([fname '.wnd']);
count = 1;
line = fgetl(fid);

while ischar(line) 
    
     if strcmp(line(1),'!')
         line = fgetl(fid);
     else
         VWT (count,1:8) = str2num(line);
         line = fgetl(fid);
         count = count +1;
     end
     
end

 fclose(fid) ;
 
 