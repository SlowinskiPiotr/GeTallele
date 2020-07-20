%% uses modified version of ternaxes and ternlable https://uk.mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot
% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot


gp=[0.05 0.05];
m1=[0.05 0.05];
m2=[0.05 0.025];

figure(1)
for dataset=22
    dataset,
    clf
    prp_vec_Tex=all_prp(dataset).prp_vec_Tex;
    ev_vec_Tex=all_prp(dataset).ev_vec_Tex;
    pop_vec_Tex=all_prp(dataset).pop_vec_Tex;
    
    y_Tex = prp_vec_Tex(:,1)*sin(deg2rad(60));
    x_Tex = prp_vec_Tex(:,2) + y_Tex*cot(deg2rad(60));
    
    clor=lines(7);
    clor=clor([3 6 5],:);
    p_vec=2;
    e_vec=2;
    subtightplot(2,3,3,gp,m1,m2)
    idx=(ev_vec_Tex==e_vec & pop_vec_Tex==p_vec);    
    
    ga_ternaxes(10,clor(1,:),clor(2,:),clor(3,:),8);
    hold on,
    if p_vec==4
        ga_ternlabel('T2+T3','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.15)
    else
        ga_ternlabel('T2','T1','N', clor(2,:),clor(1,:),clor(3,:),10,0.15)
    end
    
    if sum(idx)<5 || p_vec==2
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),60,'k','filled','markerfacealpha',0.25)
    elseif sum(idx)<20 && p_vec>2
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),40,'k','filled','markerfacealpha',0.25)
    else
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),10,'k','filled','markerfacealpha',0.25)
    end

    p_vec=3;
    e_vec=1;
    subtightplot(2,3,6,gp,m1,m2)
    idx=(ev_vec_Tex==e_vec & pop_vec_Tex==p_vec);    
    
    ga_ternaxes(10,clor(1,:),clor(2,:),clor(3,:),8);
    hold on,
    if p_vec==4
        ga_ternlabel('T2+T3','T1','N', clor(2,:),clor(1,:),clor(3,:),16,0.15)
    else
        ga_ternlabel('T2','T1','N', clor(2,:),clor(1,:),clor(3,:),16,0.15)
    end
    
    if sum(idx)<5 || p_vec==2
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),60,'k','filled','markerfacealpha',0.25)
    elseif sum(idx)<20 && p_vec>2
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),40,'k','filled','markerfacealpha',0.25)
    else
        xy=unique([x_Tex(idx),y_Tex(idx)],'rows');
        scatter(xy(:,1),xy(:,2),10,'k','filled','markerfacealpha',0.25)
    end
    
end
%
subtightplot(2,3,2,gp,m1,m2)
cla
cc=dictionaries(2,2).dict(:,:,160-99)'; 
cc(1,1)=0;
imagesc(triu(cc),[0 1])

axis square
evA_lbl={'0','A','AA'};
evB_lbl={'0','B','BB'};

set(gca,'TickDir','out','XTick',1:6,'YTick',1:6,'XTickLabel',evA_lbl,'YTickLabel',evB_lbl,'xaxisLocation','top')
xtickangle(-45)

cl=brewermap(20,'Greys');
colormap(cl)
grid on
for i1=1:3
    for i2=1:3
        if cc(i1,i2)>0
        text(i1,i2,num2str(cc(i1,i2),3),'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end
end

subtightplot(2,3,5,gp,m1,m2)
cla
cc=dictionaries(3,1).dict(:,:,3350)';
cc(1:2,1:2)=0;

imagesc(triu(cc),[0 1])

axis square
evA_lbl={'T1:0, T2:0','T1:0, T2:A','T1:A, T2:A'};
evB_lbl={'T1:0, T2:0','T1:0, T2:B','T1:B, T2:B'};

set(gca,'TickDir','out','XTick',1:6,'YTick',1:6,'XTickLabel',evA_lbl,'YTickLabel',evB_lbl,'xaxisLocation','top')
xtickangle(-45)
colormap(cl)
grid on
for i1=1:3
    for i2=1:3
        if cc(i1,i2)>0
        text(i1,i2,num2str(cc(i1,i2),3),'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end
end