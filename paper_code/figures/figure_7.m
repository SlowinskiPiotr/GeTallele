%% uses subtightplot function from https://uk.mathworks.com/matlabcentral/fileexchange/39664-subtightplot
dataset=14;
data=participant_BRCA(dataset);
all_data_pos=all_data(dataset).data(:,1:2);

chrm=1;
data_pos=all_data_pos(all_data_pos(:,1)==chrm,2);

vprs_based_on_Tex=data.all_vprs_mat_Tex;
vprs_based_on_Tex(vprs_based_on_Tex(:,10)>1e-5 & vprs_based_on_Tex(:,6)>=0.58,6)=NaN;
vprs_based_on_Tex(vprs_based_on_Tex(:,11)>1e-5 & vprs_based_on_Tex(:,7)>=0.58,7)=NaN;

vprs_based_on_Ttr=data.all_vprs_mat_Ttr;
vprs_based_on_Ttr(vprs_based_on_Ttr(:,10)>1e-5 & vprs_based_on_Ttr(:,6)>=0.58,6)=NaN;
vprs_based_on_Ttr(vprs_based_on_Ttr(:,11)>1e-5 & vprs_based_on_Ttr(:,7)>=0.58,7)=NaN;

vprs_based_on_Tex=vprs_based_on_Tex(vprs_based_on_Tex(:,1)==chrm,:);
vprs_based_on_Ttr=vprs_based_on_Ttr(vprs_based_on_Ttr(:,1)==chrm,:);

[vprs_pos_on_Tex,ims_on_Tex]=unique([vprs_based_on_Tex(:,2); vprs_based_on_Tex(:,3)+eps]);
all_vprs_Ttr_on_Tex=[vprs_based_on_Tex(:,7); vprs_based_on_Tex(:,7)];
all_vprs_Ttr_on_Tex=all_vprs_Ttr_on_Tex(ims_on_Tex);

all_vprs_Tex_on_Tex=[vprs_based_on_Tex(:,6); vprs_based_on_Tex(:,6)];
all_vprs_Tex_on_Tex=all_vprs_Tex_on_Tex(ims_on_Tex);

[vprs_pos_on_Ttr,ims_on_Ttr]=unique([vprs_based_on_Ttr(:,2); vprs_based_on_Ttr(:,3)+eps]);
all_vprs_Ttr_on_Ttr=[vprs_based_on_Ttr(:,7); vprs_based_on_Ttr(:,7)];
all_vprs_Ttr_on_Ttr=all_vprs_Ttr_on_Ttr(ims_on_Ttr);

all_vprs_Tex_on_Ttr=[vprs_based_on_Ttr(:,6); vprs_based_on_Ttr(:,6)];
all_vprs_Tex_on_Ttr=all_vprs_Tex_on_Ttr(ims_on_Ttr);

interp_vprs_Tex_on_Tex=interp1(vprs_pos_on_Tex,all_vprs_Tex_on_Tex,data_pos,'nearest','extrap');
interp_vprs_Ttr_on_Tex=interp1(vprs_pos_on_Tex,all_vprs_Ttr_on_Tex,data_pos,'nearest','extrap');
interp_vprs_Tex_on_Ttr=interp1(vprs_pos_on_Ttr,all_vprs_Tex_on_Ttr,data_pos,'nearest','extrap');
interp_vprs_Ttr_on_Ttr=interp1(vprs_pos_on_Ttr,all_vprs_Ttr_on_Ttr,data_pos,'nearest','extrap');

colr=lines(7);
gap=[0.1 0.05];
marg_h=[0.05 0.05];
marg_w=[0.075 0.025];

subtightplot(6,1,[1 2],gap,marg_h,marg_w)
cla
plot(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,1),abs(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,5)-0.5)+0.5,...
    '.','markersize',3,'color',colr(3,:))
hold on,
plot(data_pos,interp_vprs_Ttr_on_Tex,'o-','markersize',3,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
plot(data_pos,interp_vprs_Ttr_on_Ttr,'x-','markersize',3,'color',colr(3,:))
set(gca,'xTickLabel','','yTick',[0.5 0.75 1])
ylabel('VAF')
ylim([0.5 1])
legend({'VAF_{TTR}','v_{PR,TTR} in windows based on VAF_{TEX}','v_{PR,TTR} in windows based on VAF_{TTR}'})
grid on

subtightplot(6,1,3,gap,marg_h,marg_w)
cla
[~,ia] = unique(data_pos);
bar(data_pos(ia),abs(interp_vprs_Ttr_on_Tex(ia)-interp_vprs_Ttr_on_Ttr(ia)),100000,'k')
ylim([0 0.3])
set(gca,'xTickLabel','','yTick',[0 0.15 0.3])
title('Absolute difference between v_{PR,TTR}''s in windows based on VAF_{TEX} and VAF_{TTR} signals.')
grid on
ylabel('MAE')

subtightplot(6,1,[4 5],gap,marg_h,marg_w)
cla
plot(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,1),abs(participant_BRCA(dataset).chrm(chrm).data_and_stats_Tex.segment(:,4)-0.5)+0.5,...
    '.','markersize',3,'color',colr(2,:),'MarkerFaceColor',colr(2,:))
hold on,
plot(data_pos,interp_vprs_Tex_on_Tex,'o-','markersize',3,'color',colr(2,:),'MarkerFaceColor',colr(2,:))
plot(data_pos,interp_vprs_Tex_on_Ttr,'x-','markersize',3,'color',[0.5 0.5 0.5])
ylim([0.5 1])
set(gca,'xTickLabel','','yTick',[0.5 0.75 1])
legend({'VAF_{TEX}','v_{PR,TEX} in windows based on VAF_{TEX}','v_{PR,TEX} in windows based on VAF_{TTR}'})
grid on
ylabel('VAF')

subtightplot(6,1,6,gap,marg_h,marg_w)
cla
bar(data_pos(ia),abs(interp_vprs_Tex_on_Tex(ia)-interp_vprs_Tex_on_Ttr(ia)),100000,'k')
ylim([0 0.3])
set(gca,'xTickLabel','','yTick',[0 0.15 0.3])
title('Absolute difference between v_{PR,TEX}''s in windows based on VAF_{TEX} and VAF_{TTR} signals.')
xlabel('BP along a chromosome')
hold off
grid on
ylabel('MAE')

