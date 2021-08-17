function [stats]=LME_analysis(data)

[stats]=LME_prep(data);
%[lfv,stats]=parsvd(tmp.X);
%[stats]=parsig(stats);
[lfv]=mtm_svd_lfv(stats.X,2,3,1/12,0,0);
[fr, evalper]=mtm_svd_conf(stats.X,2,3,1/12,100,0,0,0,0);

stats.lfv=lfv.spectrum;
stats.f=lfv.specdomain;
stats.clv=evalper';
