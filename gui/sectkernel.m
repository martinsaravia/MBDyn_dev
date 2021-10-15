clc;   clear -handles

disp ('---------- LAUNCHING PRE SCRIPT ----------')

FE.job =  get(handles.fileedit,'String');

FE.ext =  handles.exte;

FE.rundir = handles.dire;

filename=[FE.job FE.ext]; 

secflag = zeros ()

time=cputime; disp (['Processing: ' filename])
 
[NS,ES,D,DM,LS,MS,wf,info] = sectgen(filename, 1, 1,'outer','NO', 'LINEAR');
 
plotd (handles,NS,ES,wf);
 
disp (['Sectional Properties solved in: ' num2str(cputime-time) ' seconds'])