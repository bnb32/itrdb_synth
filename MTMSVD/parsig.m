function mtmstat=parsig(mtmstat,outname);
%clvstat=sigmtmsvd1(mtmstat);
%
%Written 02/06
%
%Last Modified: 7.11.08
%
%Bootstrap significance test for results of MTM-SVD. This program uses the
%bootstrap method described by Mann and Park (1999) to estimate confidence
%intervals on the lfv spectrum.
%
%The default number of realizations is 100 (needed to compute the 99%
%confidence interval), so it runs fairly slowly. Once it has been run, the
%intervals are saved in conf.mat

X=mtmstat.X;
lfv=mtmstat.lfv;
[N,M] = size(X);
nr = 100;
t = [1:N];
nfboot = N;
ufboot = round(nfboot/2);
f0bwl = mtmstat.f0bwl;
[empty, ind] = bootstrp(nr, [], t);

disp(['Performing ' num2str(nr) ' bootstrap realizations: ']);

countDir='/home/bnb32/projects/itrdb_synth/MTMSVD/countDir';
psimple_progress(countDir,0,nr);

if isempty(gcp('nocreate'));
    ncores=feature('numcores');
    parpool(ncores);
end

parfor k = 1:nr;
    psimple_progress(countDir,k,nr)
    Zboot=X(ind(:,k),:);
    %lfvboot = mtmsvdf1(Zboot,nfboot);
    [lfvboot,~] = parsvd(Zboot,nfboot);
    lfv_set(k,:) = (lfvboot);
    f0_set(k,:) = (lfvboot(1));
     
end

uf = mtmstat.uf;
f = mtmstat.f;
ord_lfv = sort(lfv_set,'descend');
ord_f0 = sort(f0_set, 'descend');
onevec = ones(1,uf);
med_lfv = median(ord_lfv(:,2:end)')';
med_f0 = (ord_f0)';
lfvall = lfv;
conf_f = uf;

clv=confonly(mtmstat,lfvall,nr,med_lfv,med_f0,lfv_set,outname);

%clvstat.uf =uf;
%clvstat.f =f;
%clvstat.ord_lfv =ord_lfv;
%clvstat.ord_f0 =ord_f0;
%clvstat.onevec = onevec;
%clvstat.med_lfv =med_lfv;
%clvstat.med_f0=med_f0;
%clvstat.lfvall =lfvall; 
%clvstat.conf_f = conf_f;
%clvstat.clv=clv;
%clvstat.nr=nr;
%clvstat.lfv_set=lfv_set;

mtmstat.uf =uf;
mtmstat.f =f;
mtmstat.ord_lfv =ord_lfv;
mtmstat.ord_f0 =ord_f0;
mtmstat.onevec = onevec;
mtmstat.med_lfv =med_lfv;
mtmstat.med_f0=med_f0;
mtmstat.lfvall =lfvall; 
mtmstat.conf_f = conf_f;
mtmstat.clv=clv;
mtmstat.nr=nr;
mtmstat.lfv_set=lfv_set;
idx=mtmstat.f0bwl+1;
mtmstat.pvals=[clv(idx,1) clv(idx,2) clv(idx,3) clv(idx,4) clv(idx,5)];
