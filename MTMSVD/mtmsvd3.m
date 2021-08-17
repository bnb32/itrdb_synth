clear
close all

load yr_opt.mat
q1 = 1:size(V2.tsm,2);
t1 = 291:391;

q2 = find(V2.tvalid <= 1890);
t2 = 291:391;

q3 = find(V2.tvalid <= 1850);
t3 = 251:391;

q4 = find(V2.tvalid <= 1800);
t4 = 191:391;

q5 = find(V2.tvalid <= 1750);
t5 = 150:391;

q6 = find(V2.tvalid <= 1700);
t6 = 100:391;

tvalid = t4;
  qtsm = q4;

X = V2.tsm(tvalid,qtsm);
Xall = V2.tsm(tvalid,:);
time = V2.time(tvalid);
% load manndata
% X = data.tsm;
% time = 1:N;

[N,M] = size(X);
nf = 1600;

[lfv,L,U,V,mtmstat] = mtmsvdf1(X,nf);
lfvall = mtmsvdf1(Xall,nf);

f = mtmstat.f;
uf = mtmstat.uf;
plot(lfvall,':','Color',[.6 .6 .6])
hold on
plot(lfv(1:round(uf)),'k')
yrttl = [num2str(time(1)) '-' num2str(time(end))];