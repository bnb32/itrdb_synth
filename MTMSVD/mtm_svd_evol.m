function [eof]=mtm_svd_evol(sh,nw,k,dt,freq,imode,npad,lan,vw,win,inc)

win_step=round(win/dt);
inc_step=round(inc/dt);
[N,M]=size(sh);

disp('')
disp(['EVOLUTIVE RECONSTRUCTION'])
disp([' Window Size : ' num2str(win) ' years'])
disp([' Window Size : ' num2str(win_step) ' time steps'])
disp([' Overlap Size : ' num2str(inc_step) ' time steps'])
disp('')

if le(win_step,N);
    temp=sh(1:win_step-1,:);
else 
    temp=sh(1:end,:);
end    
[R,RC,RP,U,S,V]=mtm_svd_bandrecon(temp,nw,k,dt,freq,imode,npad,lan,vw);
eof=RP.eof;

%eof=mean(eof,2);

i=inc_step;
num_steps=1;
while(lt(i+win_step-1,N));
    temp=sh(i:i+win_step-1,:);
    [R,RC,RP,U,S,V]=mtm_svd_bandrecon(temp,nw,k,dt,freq,imode,npad,lan,vw);
    eof=eof+RP.eof;
    i=i+inc_step;
    num_steps=num_steps+1;
end

eof=eof/num_steps;
