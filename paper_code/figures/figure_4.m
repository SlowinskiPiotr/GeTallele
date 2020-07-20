%%
% first run find_confidence_intervals_for_vpr_estimation_from_VAF.m
%%
% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
figure(1)
gp=[0.075 0.075];
mbt=[0.075 0.05];
mlr=[0.075 0.025];

subtightplot(2,3,1,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode4(i,:),[2.5 25 50 75 97.5]);
    h(1)=plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1);
    hold on
    h(2)=plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2);
    h(3)=plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)
title('Based on 4 VAFs')
ylabel('Estimated value')
legend(h,'95%','50%','med','Location','northwest')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
drawnow
hold off

subtightplot(2,3,2,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode10(i,:),[2.5 25 50 75 97.5]);
    plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1)
    hold on
    plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2)
    plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)

title('Based on 10 VAFs')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
drawnow
hold off

subtightplot(2,3,3,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode20(i,:),[2.5 25 50 75 97.5]);
    plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1)
    hold on
    plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2)
    plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)

title('Based on 20 VAFs')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
drawnow
hold off

subtightplot(2,3,4,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode50(i,:),[2.5 25 50 75 97.5]);
    plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1)
    hold on
    plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2)
    plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)

title('Based on 50 VAFs')
xlabel('True value')
ylabel('Estimated value')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
drawnow
hold off

subtightplot(2,3,5,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode150(i,:),[2.5 25 50 75 97.5]);
    plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1)
    hold on
    plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2)
    plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)

title('Based on 150 VAFs')
xlabel('True value')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
drawnow
hold off

subtightplot(2,3,6,gp,mbt,mlr)
cla
for i=50:99
    prcts=prctile(mode300(i,:),[2.5 25 50 75 97.5]);
    plot((i/100)*ones(1,2),prcts([1 5]),'color',[0.5 0.5 0.5],'linewidth',1)
    hold on
    plot((i/100)*ones(1,2),prcts([2 4]),'color',[0.25 0.25 0.25],'linewidth',2)
    plot((i/100),prcts(3),'ro','markersize',2,'markerfacecolor','r');
    axis square
    grid on,
end
plot([0.5 1],[0.5 1],'k','linewidth',0.5)

title('Based on 300 VAFs')
xlabel('True value')
set(gca,'TickDir','out','xtick',0.5:0.1:1,'ytick',0.5:0.1:1,'MinorGridLineStyle','-')
ax=gca;
ax.XAxis.MinorTickValues=0.55:0.1:0.95;
ax.YAxis.MinorTickValues=0.55:0.1:0.95;
grid minor
axis([0.5 1 0.5 1])
hold off