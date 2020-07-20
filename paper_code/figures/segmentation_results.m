out_Tex=[];
for dataset=1:72
    out_loop=participant_BRCA(dataset).all_vprs_mat_Tex;
    out_loop(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out_loop,1),1);
    out_loop(:,26)=dataset;
    out_Tex=[out_Tex; out_loop];
end
out_Tex=out_Tex(:,[26 1:25]);
out_Tex=sortrows(out_Tex,[1 2 3]);

out_analysis=out_Tex;
clc

disp('All windows:')
size(out_analysis,1)

idx_non=out_analysis(:,12)>1e-5 & out_analysis(:,8)>=0.58 |...
        out_analysis(:,11)>1e-5 & out_analysis(:,7)>=0.58;

disp('Excluded windows:')
sum(idx_non)

disp('% of bp in excluded windows:')
sum(out_analysis(idx_non,4)-out_analysis(idx_non,3))/sum(out_analysis(:,4)-out_analysis(:,3))

disp('% of data points in excluded windows:')
sum(out_analysis(idx_non,20))/sum(out_analysis(:,20))

out_analysis(idx_non,:)=[]; %remove 294 windows

disp('% of statistically concordant windows')
idx=out_analysis(:,18)>1e-5;
s1=sum(idx)/size(out_analysis,1)

disp('same vPR but statistically different')
idx=out_analysis(:,7)==out_analysis(:,8) & out_analysis(:,18)<=1e-5;
sum(idx),

disp('# of discordant and statistically different')
idx=out_analysis(:,7)~=out_analysis(:,8) & out_analysis(:,18)<=1e-5; 
sum(idx)
disp('% of discordant and statistically different windows')
sum(idx)/size(out_analysis,1)

disp('v_{PR,TTR}>v_{PR,TEX}')
idx=(out_analysis(:,7)>out_analysis(:,8)) & out_analysis(:,7)~=out_analysis(:,8) & out_analysis(:,18)<=1e-5;
sum(idx)