function LME_plot(stats)

outname='./results/LME';

plot_signals(stats,outname,1850,2005,3,0,0.5,1/12);

clv=stats.clv;

fig=figure(1);
plot(stats.f,stats.lfv,'k','linewidth',1.5);
hold on
%text(.505,p999(end),'99.9%','Fontsize',8)
text(.505,clv(end,1),'99%','Fontsize',12)
text(.505,clv(end,2),'95%','Fontsize',12)
text(.505,clv(end,3),'90%','Fontsize',12)
text(.505,clv(end,4),'80%','Fontsize',12)
text(.505,clv(end,5),'Median','Fontsize',12)
plot(stats.f,[clv(:,1),clv(:,2),clv(:,3),clv(:,4),clv(:,5)],'k:');
lfv_plot=[outname '_lfv.png'];
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
axis([0 0.5 -inf inf]);
saveas(fig,lfv_plot);
hold off
clf,cla
disp(['saved lfv spectrum as ' lfv_plot]) 
