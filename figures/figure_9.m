%% uses function distinguishable_colors from https://uk.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
% uses function Violin from https://github.com/bastibe/Violinplot-Matlab
% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

figure(1)
gp=[0.075 0.05];
m1=[0.065 0.01];
m2=[0.065 0.01];
VBP_out=NaN(1,72);
VBP_median_out=NaN(1,72);

for i=1:72
    purity_estim(i,:)=participant_BRCA(i).purity;
end

clf
dcolors = distinguishable_colors(12,[0.75 0.75 0.75]);
dcolors=dcolors(1:2:12,:);
%A
subtightplot(3,3,1:3,gp,m1,m2)
cla
for i=1:72
    i
    if ~isempty(VBP_all(i).pur_dist)
        vbp_data_out=[];
        for p_i=2:4
            for e_i=1:4
                if ~isempty(VBP_all(i).pur_dist(p_i,e_i).pr)
                    vbp_data=VBP_all(i).pur_dist(p_i,e_i).pr;
                    vbp_data_out=vbp_data;
                    VBP_median_out(i)=median(vbp_data);
                    VBP_out(i)=min(vbp_data);
                    break
                end
            end
        end


    vl=Violin(vbp_data_out,i,'ShowData',false,'ViolinColor',dcolors(6,:),'Width',0.4);
    set(vl.MedianPlot,'HandleVisibility','off')
    set(vl.BoxPlot,'HandleVisibility','off')
    
    if ~isstruct(vl.WhiskerPlot)
        set(vl.WhiskerPlot,'HandleVisibility','off')
    end
    set(vl.ScatterPlot,'HandleVisibility','off')
    set(vl.ViolinPlot,'HandleVisibility','off')
    set(vl.NotchPlots,'HandleVisibility','off')
    set(vl.MeanPlot,'HandleVisibility','off')

    hold on
    end
end
set(gca,'TickDir','out','Xtick',1:72,'XTickLabel',[])
xlim([0 73])
grid on,
 
symbs={'o','o','o','o','o'};
title_str={'ESTIMATE','ABSOLUTE','LUMP','IHC','CPE'};
plot(1:72,VBP_median_out,'o','MarkerFaceColor','w','MarkerEdgeColor','k','markersize',2)

for i=1:5
plot(1:72,purity_estim(:,i),symbs{i},'Color',dcolors(i,:),'MarkerFaceColor',dcolors(i,:),'markersize',3,'HandleVisibility','off')
end

plot(1:72,VBP_out,'x','color',[0.4 0.4 0.4],'markersize',4,'linewidth',1)
legend({'median','VBP estimate'},'Location','best','box','off')
set(gca,'color',[0.85 0.85 0.85])
xlabel('sample')
ylabel('purity')
hold off

%B1-B5
for i=1:5
    subtightplot(3,3,i+3,gp,m1,m2)
    cla
    plot(VBP_out,purity_estim(:,i),symbs{i},'color',dcolors(i,:),'markersize',3,'MarkerFaceColor',dcolors(i,:));
    lsline
    [r,p]=corr(VBP_out',purity_estim(:,i),'rows','complete','type','Pearson');
    hold on
    plot([0 1],[0 1],'-','Color',[.5 .5 .5])
    title([title_str{i} ': r=' num2str(r,2) '; p<' num2str(p,2)],'FontWeight','normal')
    axis([0 1 0 1])
    axis square
        set(gca,'TickDir','out','XTick',0:0.5:1,'YTick',0:0.5:1)

    if i==1 || i==4
        ylabel('purity')
    else
        set(gca,'YTickLabel',[])
    end
    if i<3 
        set(gca,'XTickLabel',[])
    else
        xlabel('V_{PR} based purity (VBP)')
    end
    grid on
    box off
    hold  off
end

%C
[cc,p]=corr([VBP_out' purity_estim],'rows','pairwise','type','Pearson');
cc(p>0.05)=0;

subtightplot(3,3,9,gp,m1,m2)
imagesc(triu(cc)',[0.3 1])
axis square
method_lbl={'VBP','EST','ABS','LUMP','IHC','CPE'};
set(gca,'TickDir','out','XTick',1:6,'YTick',1:6,'XTickLabel',method_lbl,'YTickLabel',method_lbl)
colormap(brewermap(15,'Greys'))
grid on
for i1=1:6
    for i2=i1+1:6
        if cc(i1,i2)>0
        text(i1,i2,num2str(cc(i1,i2),2),'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end
end