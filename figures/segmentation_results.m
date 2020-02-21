%%
out_Tex=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_modes_mat_Tex;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Tex=[out_Tex; out];
end
out_Tex=out_Tex(:,[26 1:25]);
out_Tex=sortrows(out_Tex,[1 2 3]);
%csvwrite('all_modes.csv',out)
%%
out_Ttr=[];
for dataset=1:72
    out=participant_BRCA(dataset).all_modes_mat_Ttr;
    out(:,21:25)=repmat(participant_BRCA(dataset).purity,size(out,1),1);
    out(:,26)=dataset;
    out_Ttr=[out_Ttr; out];
end
out_Ttr=out_Ttr(:,[26 1:25]);
out_Ttr=sortrows(out_Ttr,[1 2 3]);
%csvwrite('all_modes.csv',out)
%%
out=out_Tex;
clc

disp('All windows:')
size(out,1)
idx_non=out(:,12)>1e-5 & out(:,8)>=0.58 |...
         out(:,11)>1e-5 & out(:,7)>=0.58;

disp('Excluded windows:')
sum(idx_non)

disp('% of bp in excluded windows:')
sum(out(idx_non,4)-out(idx_non,3))/sum(out(:,4)-out(:,3))
%sum(out(~idx_non,4)-out(~idx_non,3))/sum(out(:,4)-out(:,3))

disp('% of data points in excluded windows:')
sum(out(idx_non,20))/sum(out(:,20))
%sum(out(~idx_non,20))/sum(out(:,20))

out(idx_non,:)=[];

idx05_Tex=out(:,11)>1e-5;
out(idx05_Tex,7)=0.5;
idx05_Tex=out(:,12)>1e-5;
out(idx05_Tex,8)=0.5;

size(out,1)
disp('out:')
disp('concordant same and one=0.5')
idx=(out(:,7)==0.5 | out(:,8)==0.5) & out(:,18)>1e-5; % concordant =0.5
s1=sum(idx)/size(out,1)

disp('concordant both>0.5')
idx=(out(:,7)>0.5 & out(:,8)>0.5) & out(:,18)>1e-5; % concordant >0.5
s2=sum(idx)/size(out,1)

disp('concordant but statistically different')
idx=out(:,7)==out(:,8) & out(:,18)<=1e-5; % concordant but statistically different
s3=sum(idx)/size(out,1)

disp('% of concordant windows')
s1+s2+s3

disp('EXTRA concordant both<0.6')
idx=out(:,7)<0.6 & out(:,8)<0.6 & out(:,18)>1e-5; % concordant >0.5
sum(idx)/size(out,1)
disp('EXTRA concordant and identical')
idx=out(:,7)==out(:,8) & out(:,18)>1e-5; % concordant and identical
sum(idx)/size(out,1)
disp('EXTRA concordant both=0.5')
idx=(out(:,7)==0.5 & out(:,8)==0.5) & out(:,18)>1e-5; % concordant =0.5
sum(idx)/size(out,1)

disp('discordant and statistically different')
idx=out(:,7)~=out(:,8) & out(:,18)<=1e-5; % discordant and statistically different
sum(idx)
disp('% of discordant windows')
sum(idx)/size(out,1)

disp('EXTRA discordant and statistically different Tex>0.5')
idx=out(:,7)>0.5 & out(:,7)~=out(:,8) & out(:,18)<=1e-5;
sum(idx)/size(out,1)
disp('EXTRA discordant and statistically different Tex=0.5')
idx=out(:,7)==0.5 & out(:,7)~=out(:,8) & out(:,18)<=1e-5;
sum(idx)/size(out,1)
disp('EXTRA discordant and statistically different Ttr=0.5')
idx=out(:,7)>0.5 & out(:,8)==0.5  & out(:,18)<=1e-5;
sum(idx)/size(out,1)

disp('EXTRA V_{prTTR}>V_{prTEX}')
idx=(out(:,7)>out(:,8)) & out(:,7)~=out(:,8) & out(:,18)<=1e-5;
sum(idx)