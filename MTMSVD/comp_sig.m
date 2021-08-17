function tdiff=comp_sig(file1,file2,range1,range2)

outname=[file1(1:end-20) '_' range1 'vs' range2 '_compsig.png'];

s1=load(file1);
s2=load(file2);

l1=s1.stats.lfv;
f1=s1.stats.f(1:length(l1));

l2=s2.stats.lfv;
f2=s2.stats.f(1:length(l2));

if length(l1)>length(l2)
    l2=interp1(f2,l2,f1,'linear','extrap');
    f2=f1;
else
    l1=interp1(f1,l1,f2,'linear','extrap');
    f1=f2;
end   

tdiff=0;

for i=1:length(l1)
    diff=abs(l1(i)-l2(i));
    avg=(l1(i)+l2(i))/2.0;
    if ~isnan(diff) && ~isnan(avg)
        tdiff=tdiff+(diff/avg);
    end
end

tdiff=tdiff/length(l1);

disp(['diff= ',num2str(tdiff)])

fig=figure(3);
cla

plot(f1,l1,'b', 'linewidth', 1.5)
hold on
plot(f2,l2,'r', 'linewidth', 1.5)

ax=gca;
ax.FontSize=18;
legend(range1, range2)
xl=get(gca,'XTickLabel');
yl=get(gca,'YTickLabel');
set(gca,'XTickLabel',xl,'FontSize',14,'FontWeight','normal');
set(gca,'YTickLabel',yl,'FontSize',14,'FontWeight','normal');
xlabel('Frequency (Cycles/year)','FontSize',18,'FontWeight','Bold');
ylabel('Amplitude','FontSize',18,'FontWeight','Bold');
set(gca,'XTickLabelMode','auto');
set(gca,'YTickLabelMode','auto');
title('Local Fractional Variance');
saveas(fig,outname)

end
