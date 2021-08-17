function ts_plot_gen(ts,outfile,title_name,syear,eyear);

fig=figure(1);
cla,clf;

ts=mean(ts,1);

[M,N]=size(ts);

times=linspace(syear,eyear,length(ts));

plot(times,ts,'linewidth',1.5);

ax=gca;
ax.FontSize=16;
%plot(amp(1:end));
%plot(exps(1:end));
%plot(w(:,1));
%plot(w(:,2));
%plot(w(:,3));
xlabel('Year','FontSize',14,'FontWeight','Bold');
ylabel('Amplitude','FontSize',14,'FontWeight','Bold');
title(title_name);
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

axis([-inf inf -inf inf]);
saveas(fig,outfile);
disp(['saved time series plot as ' outfile])

