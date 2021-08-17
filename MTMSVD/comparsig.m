% load clv_q1t1;
% a = [p99(1) p99(end)];
% load clv_q1t4
% b = [p99(1) p99(end)];
% load clv_q4t4
% c = [p99(1) p99(end)];
% load clv_q5t5
% d = [p99(1) p99(end)];
% load clv_q6t6
% e = [p99(1) p99(end)];
% 
% plot([a' b' c' d' e'],'.','LineWidth',2)
% legend('a','b','c','d','e')

clear
close all

load yr_opt.mat
q1 = 1:size(V2.tsm,2);
t1 = 291:391;

q2 = find(V2.tvalid <= 1890);
t2 = 291:391;

q3 = find(V2.tvalid <= 1850);
t3 = 250:391;

q4 = find(V2.tvalid <= 1800);
t4 = 191:391;

q5 = find(V2.tvalid <= 1750);
t5 = 150:391;

q6 = find(V2.tvalid <= 1700);
t6 = 100:391;

tvalid = t4;
  qtsm = q4;

X = V2.tsm(tvalid,qtsm);
time = V2.time(tvalid);
% load manndata
% X = data.tsm;
% time = 1:N;

[N,M] = size(X);
nf = 4*N;

[lfv,L,U,V,mtmstat] = mtmsvdf1(X,nf);
f = mtmstat.f;
uf = mtmstat.uf;
plot(lfv(1:round(uf)))