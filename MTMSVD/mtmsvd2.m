clear
close all

load yr_opt.mat
tvalid = 1:400;
X = V2.tsm(tvalid,:);
time = V2.time(tvalid);

% load manndata
% X = data.tsm;
% time = 1:N;

[N,M] = size(X);
nf = N*4;

[lfv,L,U,V,mtmstat] = mtmsvdf1(X);

f = mtmstat.f;
uf = mtmstat.uf;
plot(f(1:round(uf)),lfv(1:round(uf)))