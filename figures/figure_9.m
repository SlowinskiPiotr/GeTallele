%% uses function distinguishable_colors from https://uk.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
% uses function Violin from https://github.com/bastibe/Violinplot-Matlab
figure
title_str={'ESTIMATE','ABSOLUTE','LUMP','IHC','CPE'};
for i=1:5
    subtightplot(2,3,i,[0.025 0.05],[0.065 0.01],[0.065 0.01])
    plot(VBP_out',purity_estim(:,i),'x','markersize',3,'MarkerFaceColor','r');%lines(1))
    lsline
    [r,p]=corr(VBP_out',purity_estim(:,i),'rows','complete','type','Pearson');
    hold on
    plot([0 1],[0 1],'-','Color',[.5 .5 .5])
    title([title_str{i} ': r=' num2str(r,2) '; p<' num2str(p,2)])
    axis([0 1 0 1])
    axis square
        set(gca,'TickDir','out')

    if i==1 || i==4
        ylabel('purity')
    else
        set(gca,'YTickLabel',[])
    end
    if i<3 
        set(gca,'XTickLabel',[])
    else
        xlabel('V_{PR} based purity')
    end
    box off
    hold  off
end

[cc,p]=corr([VBP_out' purity_estim],'rows','pairwise','type','Pearson');

triu(cc)
triu(p)

%%
figure(1)
gp=[0.075 0.05];
m1=[0.065 0.01];
m2=[0.065 0.01];
VBP_out=NaN(1,72);
VBP_out_median=NaN(1,72);

clf
dcolors = distinguishable_colors(12,[0.75 0.75 0.75]);
dcolors=dcolors(1:2:12,:);
subtightplot(3,3,[1:3],gp,m1,m2)
cla
for i=1:72
    i
    if ~isempty(VBP_all(i).pur_dist)%sum(all_pr_idx==i)>0
        vbp_data_out=1;
        b=0;
        for e_i=1:4
            for p_i=2:4
                if ~isempty(VBP_all(i).pur_dist(p_i,e_i).pr)
                    vbp_data=VBP_all(i).pur_dist(p_i,e_i).pr;
                    vbp_data_out=vbp_data;
                    VBP_out_median(i)=median(vbp_data);
                    VBP_out(i)=min(vbp_data);
                    b=1;
                    break
%                     if min(vbp_data)<min(vbp_data_out)
%                        VBP_out(i)=min(vbp_data);
%                        VBP_out_median(i)=median(vbp_data);
%                        vbp_data_out=vbp_data;
%                     end
                end
            end
            if b==1
                b=0;
                break
            end
        end
    %DBP_out(i)=min(dbp_data);
    %DBP_out_median(i)=median(dbp_data);

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

%     vl.BoxColor='none';
%     vl.MedianPlot.Marker='none';
%     vl.EdgeColor=vl(i).ViolinColor;
    hold on
    end
end
set(gca,'TickDir','out','Xtick',1:72,'XTickLabel',[])
xlim([0 73])

grid on,
 

%symbs={'s','+','d','*','p'};
symbs={'o','o','o','o','o'};
title_str={'ESTIMATE','ABSOLUTE','LUMP','IHC','CPE'};
plot(1:72,VBP_out_median,'o','MarkerFaceColor','w','MarkerEdgeColor','k','markersize',2)

for i=1:5
plot(1:72,purity_estim(:,i),symbs{i},'Color',dcolors(i,:),'MarkerFaceColor',dcolors(i,:),'markersize',3,'HandleVisibility','off')
end
plot(1:72,VBP_out,'x','color',[0.4 0.4 0.4],'markersize',4,'linewidth',1)

legend({'median','VBP estimate'},'Location','best','box','off')
set(gca,'color',[0.85 0.85 0.85])
xlabel('sample')
ylabel('purity')
hold off
for i=1:5
    subtightplot(3,3,i+3,gp,m1,m2)
    cla
    plot(purity_estim(:,2),purity_estim(:,i),symbs{i},'color',dcolors(i,:),'markersize',3,'MarkerFaceColor',dcolors(i,:));%lines(1))
    lsline
    [r,p]=corr(purity_estim(:,2),purity_estim(:,i),'rows','complete','type','Pearson');
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
% 
% subtightplot(2,3,6,[0.025 0.05],[0.065 0.01],[0.065 0.01])
%     plot(DBP_Tex_end,DBP_Ttr_end,'x','markersize',3,'MarkerFaceColor','r');%lines(1))
% lsline
%     hold on
%     plot([0 1],[0 1],'-','Color',[.5 .5 .5])
%         axis([0 1 0 1])
%     axis square
%         set(gca,'TickDir','out')
%         box off
%         hold off
%     [r,p]=corr(DBP_Ttr_end',DBP_Tex_end','rows','complete','type','Pearson');
%     title(['DBP_{Ttr}: r=' num2str(r,2) '; p<' num2str(p,2)])

[cc,p]=corr([VBP_out_median' purity_estim],'rows','pairwise','type','Pearson');
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
    



%%
DBP_from=NaN(1,72);
for i=1:72
    i
    if sum(all_pr_idx==i)>0
    dummy=unique(all_pr(all_pr_idx==i));
    DBP_from(i)=dummy(end);
    end
end
figure
plot(DBP_from,DBP,'o')