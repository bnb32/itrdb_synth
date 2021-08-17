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
[N,M] = size(V2.tsm);
nf = 4*N;
f = (1:nf)/nf;
uf = round(nf/2);

q{1} = 1:size(V2.tsm,2);
t{1} = 291:391;

q{2} = find(V2.tvalid <= 1890);
t{2} = 291:391;

q{3} = find(V2.tvalid <= 1850);
t{3} = 250:391;

q{4} = find(V2.tvalid <= 1800);
t{4} = 191:391;

q{5} = find(V2.tvalid <= 1750);
t{5} = 150:391;

q{6} = find(V2.tvalid <= 1700);
t{6} = 100:391;

k = [1 2 1 2 1 2];
for i = 1:6
    clear X Xall lfvx lfvall
    X = V2.tsm(t{i},q{i});
    Xall = V2.tsm(t{i},:);
    
    lfvx = mtmsvdf1(X,nf);
    lfvall = mtmsvdf1(Xall,nf);
    
    figure(ceil(i/2))
    subplot(2,1,k(i))
    plot(f(1:uf), lfvall,'Color',[.7 .7 .7])
    hold on
    plot(f(1:uf), lfvx, 'k')
end