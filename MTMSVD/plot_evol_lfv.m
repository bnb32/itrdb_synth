function plot_evol_lfv(yrs,freqs,lfv,outname);

fig=figure(1);
pcolor(yrs,freqs,lfv');
shading('interp');
colorbar;
axis([-inf inf 0 inf]);
ax=gca;
ax.FontSize=16;
xlabel('Year','FontSize',14,'FontWeight','Bold');
ylabel('Frequency (cyc/yr)','FontSize',14,'FontWeight','Bold');
xl = get(gca,'XLabel');
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 14)
set(xl, 'FontSize', xlFontSize);
yl = get(gca,'XLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 14)
set(yl, 'FontSize', ylFontSize);
title('Evolutive LFV Spectrum');

saveas(fig,outname);
disp(['saved figure as ' outname]);
