function clvstat=sigmtmsvd1(X,lfv,mtmstat);
%clvstat=sigmtmsvd1(X,lfv,mtmstat);
%
%Written 02/06
%
%Last Modified: 7.11.08
%
%Bootstrap significance test for results of MTM-SVD. This program uses the
%bootstrap method described by Mann and Park (1999) to estimate confidence
%intervals on the lfv spectrum.
%
%The default number of realizations is 1000 (needed to compute the 99%
%confidence interval), so it runs fairly slowly. Once it has been run, the
%intervals are saved in conf.mat

[N,M] = size(X);
nr = 100;
t = [1:N];
nfboot = N;
ufboot = round(nfboot/2);
f0bwl = mtmstat.f0bwl;
[empty, ind] = bootstrp(nr, [], t);

h = waitbar(0, ['Performing ' num2str(nr) ' bootstrap realizations']);
for k = 1:nr;
    clear lfvboot Zboot
    Zboot=X(ind(:,k),:);
    lfvboot = mtmsvdf1(Zboot,nfboot);
    
    lfv_set(k,:) = (lfvboot);
    f0_set(k,:) = (lfvboot(1));
    
    waitbar(k/nr,h);
end
close(h);
uf = mtmstat.uf;
f = mtmstat.f;
ord_lfv = sort(lfv_set,'descend');
ord_f0 = sort(f0_set, 'descend');
onevec = ones(1,uf);
med_lfv = median(ord_lfv(:,2:end)')';
med_f0 = (ord_f0)';
lfvall = lfv;
conf_f = uf;

clv=confonly(mtmstat,lfvall,nr,med_lfv,med_f0);

clvstat.uf =uf;
clvstat.f =f;
clvstat.ord_lfv =ord_lfv;
clvstat.ord_f0 =ord_f0;
clvstat.onevec = onevec;
clvstat.med_lfv =med_lfv;
clvstat.med_f0=med_f0;
clvstat.lfvall =lfvall; 
clvstat.conf_f = conf_f;
clvstat.clv=clv;
clvstat.nr=nr;

