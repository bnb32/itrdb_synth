function [eof]=mtm_svd_evol_lfv(sh,dt,freqs,win,inc)

ws=round(win/dt);
hw=round(ws/2);
step=round(inc/dt);
[N,M]=size(sh);

disp('')
disp(['EVOLUTIVE EOF CONSTRUCTION'])
disp([' Window Size : ' num2str(win) ' years'])
disp([' Window Size : ' num2str(ws) ' time steps'])
disp('')


cn=[hw:step:N-hw];

for i=1:length(cn)
    yrs(i)=cn(i)*dt;
    data=sh(cn(i)-hw+1:cn(i)+hw,:);
    [R,RC,RP,~]=mtm_svd_bandrecon(data,2,3,dt,freqs,1,0,2,0);
    eof(i,:)=RP.eof;
end

eof=mean(eof,1);
