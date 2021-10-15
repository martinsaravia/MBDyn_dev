disp ('---------- LAUNCHING POST SCRIPT ----------')


if get(handles.savebox1,'Value') == 1
    hwresexp (X,FE,E,N,U,Vel,Acc,AV,'dyn')
    save(FE.job)
    disp ('Results Saved...')
end

