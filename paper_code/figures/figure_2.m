%% Figure 2
%  uses subtightplot from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

gap=[0.05 0.025];
m1=[0.1 0.05];
m2=[0.075 0.025];

subtightplot(2,4,[1:3 5:7],gap,m1,m2)
sample=49;

matrix=all_data(sample).data;
layer=3;
l1=layer+6;
l2=l1+4;
idx_c=find(matrix(:,1)<=22);
VAF_vector=matrix(idx_c,l1)./(matrix(idx_c,l1)+matrix(idx_c,l2));
plot(VAF_vector,'k.')
axis([1,numel(VAF_vector) -0.01 1.01])
set(gca,'ytick',0:0.1:1,'yticklabel',0:0.1:1,'XTickLabel','','tickdir','out','YGrid','on')
xlabel('samples along chromosomes (1-22)')
ylabel('VAF')
text(0,1.05,'A','Fontweight','bold')


[bin_edges,bin_centres]=farey_bins(100);

nVAF=histcounts(VAF_vector,bin_edges);
subtightplot(2,4,[4 8],gap,m1,m2)
barh(bin_centres,nVAF,25,'facecolor','k','edgecolor','none')
set(gca,'ytick',0:0.1:1,'yticklabel','','XTickLabel','','tickdir','out','YGrid','on')
xlabel('Frequency')
ylim([-0.01 1.01])
text(0,1.05,'B','Fontweight','bold')


