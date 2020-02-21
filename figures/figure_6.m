%% uses modified version of ternaxes and ternlable https://uk.mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot
% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

gp=[0.05 0.05];
m1=[0.05 0.05];
m2=[0.05 0.025];

dataset=22
clf
prp_vec_Tex=all_prp(dataset).prp_vec_Tex;
ev_vec_Tex=all_prp(dataset).ev_vec_Tex;
pop_vec_Tex=all_prp(dataset).pop_vec_Tex;

y_Tex = prp_vec_Tex(:,1)*sin(deg2rad(60));
x_Tex = prp_vec_Tex(:,2) + y_Tex*cot(deg2rad(60));

clor=lines(7);
clor=clor([3 6 5],:);

for p_vec=2:4
    for e_vec=1:4
        
        subtightplot(3,4,e_vec+4*(p_vec-2),gp,m1,m2)
        cla
        axis off
        if sum(ev_vec_Tex==e_vec & pop_vec_Tex==p_vec)>0
            
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
            text(0.5,0.965,['T_{pl}: ' num2str(p_vec) ' ev: 0-' num2str(e_vec)],...
                'FontWeight','bold','HorizontalAlignment','Right','FontSize',10)
            
        end
    end
end



