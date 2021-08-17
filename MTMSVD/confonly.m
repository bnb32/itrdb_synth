function [clv] = confonly(mtmstat,lfvall,nr,med_lfv,med_f0,lfv_set,outname)

f = mtmstat.f;
lfv = mtmstat.lfv;
nf = length(lfv);
uf = mtmstat.uf;
fuf = mtmstat.f(1:uf);
f0bwl = mtmstat.f0bwl;

p999(1:nf) = med_lfv(ceil(.001*nr));
p999(1:f0bwl) = med_f0(ceil(.001*nr));

p99(1:nf) = med_lfv(ceil(.01*nr));
p99(1:f0bwl) = med_f0(ceil(.01*nr));

p95(1:nf) = med_lfv(ceil(.05*nr));
p95(1:f0bwl) = med_f0(ceil(.05*nr));

p90(1:nf) = med_lfv(ceil(.1*nr));
p90(1:f0bwl) = med_f0(ceil(.1*nr));

p80(1:nf) = med_lfv(ceil(.2*nr));
p80(1:f0bwl) = med_f0(ceil(.2*nr));

p50(1:nf) = median(med_lfv);
p50(1:f0bwl) = median(med_f0');


fig=figure(2);
cla

med_len=length(median(lfv_set(:,1:end)));

plot(f(1:nf),lfv(1:nf),'k','linewidth',1.5)
hold on
%plot(f(1:nf),lfvall,':','Color',[.6 .6 .6])
%plot(f(1:med_len),median(lfv_set(:,1:end)),'r','linewidth',1.5)
%legend('All Data', 'Restr')

%plot(fuf, [p999' p99' p95' p90' p80' p50'],'k:')
plot(fuf, [p99' p95' p90' p80' p50'],'k:')
%text(.505,p999(end),'99.9%','Fontsize',8)
text(.505,p99(end),'99%','Fontsize',12)
text(.505,p95(end),'95%','Fontsize',12)
text(.505,p90(end),'90%','Fontsize',12)
text(.505,p80(end),'80%','Fontsize',12)
text(.505,p50(end),'Median','Fontsize',12)

ax=gca;
ax.FontSize=18;
xl=get(gca,'XTickLabel');
yl=get(gca,'YTickLabel');
set(gca,'XTickLabel',xl,'FontSize',14,'FontWeight','normal')
set(gca,'YTickLabel',yl,'FontSize',14,'FontWeight','normal')
xlabel('Frequency (Cycles/year)','FontSize',18,'FontWeight','Bold')
ylabel('Amplitude','FontSize',18,'FontWeight','Bold')
set(gca,'XTickLabelMode','auto')
set(gca,'YTickLabelMode','auto')
title('Local Fractional Variance')

axis([0 0.505 -inf inf]);

saveas(fig,outname);
disp(['saved confidence plot as ' outname])

%clv = [p999' p99' p95' p90' p80' p50'];
clv = [p99' p95' p90' p80' p50'];
