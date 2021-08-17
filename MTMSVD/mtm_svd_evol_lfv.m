function [lfv,yrs,freqs]=mtm_svd_evol_lfv(sh,dt,win)

ws=round(win/dt);
hw=round(ws/2);
[N,M]=size(sh);

disp('')
disp(['EVOLUTIVE LFV CONSTRUCTION'])
disp([' Window Size : ' num2str(win) ' years'])
disp([' Window Size : ' num2str(ws) ' time steps'])
disp('')


cn=[hw:N-hw];

for i=1:length(cn)
    yrs(i)=cn(i)*dt;
    [tmp]=mtm_svd_lfv(sh(cn(i)-hw+1:cn(i)+hw,:),2,3,dt,0,0);
    lfv(i,:)=tmp.spectrum;
end

freqs=tmp.specdomain;
