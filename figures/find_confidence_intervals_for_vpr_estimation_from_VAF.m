bin_edges=farey_bins(1000);
bin_edges=bin_edges(152097:end);

tic
[model_05,models_cdf,all_models]=gen_ideal_hist(all_data(1).data,bin_edges);
ideal_samples_Ttr=models_cdf(3).models;
toc

layer=3;
matrix=all_data(3).data;
l1=layer+6;
l2=l1+4;
idx_c=find(matrix(:,1)<=22);
counts_vector=matrix(idx_c,l1)+matrix(idx_c,l2);
nb_reads=numel(counts_vector);    

for smp_size=[2 5 10 25 75 150]

modes=[];
for k=50:99
[k, smp_size]
tic,
parfor bts=1:1000
    surr1=[];
   
    while numel(surr1)<smp_size+1
    p=k/100;
    p=p*ones(1,smp_size);
    new_count_idx=randi(nb_reads,1,smp_size);
    new_count=counts_vector(new_count_idx);
    y = binornd(new_count(:),p(:));
    sur=y(:)./new_count(:);
    surr1=[surr1; sur];
    end
    
    surr2=[];
    while numel(surr2)<smp_size+1
    p=1-k/100;
    p=p*ones(1,smp_size);
    new_count_idx=randi(nb_reads,1,smp_size);
    new_count=counts_vector(new_count_idx);
    y = binornd(new_count(:),p(:));
    sur=y(:)./new_count(:);
    surr2=[surr2; sur];
    end 

    smpl=[surr1; surr2];
    smpl=abs(smpl-0.5)+0.5;
    
    %find mode
    Nx=sparse(histc(smpl,bin_edges));
    lx=sum(Nx);
    csx=cumsum(Nx/lx)';
    csx_mat=repmat(csx,51,1);

    emd_v=sum(abs(csx_mat-ideal_samples_Ttr),2);

    [~,i_mode]=min(emd_v);
    modes(k,bts)=(i_mode+49)/100;
end
toc,
end
eval(['mode' num2str(smp_size*2) '=modes;'])
end