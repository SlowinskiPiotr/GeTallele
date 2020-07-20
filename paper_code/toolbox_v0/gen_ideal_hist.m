function [model_05,models_cdf,all_models]=gen_ideal_hist(matrix,bin_edges)

% comment

for layer=4:-1:1
    
    l1=layer+6;
    l2=l1+4;
    
    idx_c=find(matrix(:,1)<=22);
    counts_vector=matrix(idx_c,l1)+matrix(idx_c,l2);
    
    model_samples = struct('smpl', cell(1,51));
    
    for k=50:100
        nb_reads=numel(counts_vector);
        
        p=k/100;
        p=p*ones(1,5000);
        new_count_idx=randi(nb_reads,1,5000);
        new_count=counts_vector(new_count_idx);
        y = binornd(new_count(:),p(:));
        surr1=y(:)./new_count(:);
        
        p=1-k/100;
        p=p*ones(1,5000);
        new_count_idx=randi(nb_reads,1,5000);
        new_count=counts_vector(new_count_idx);
        y = binornd(new_count(:),p(:));
        surr2=y(:)./new_count(:);
        
        model_samples(k-49).smpl=[surr1; surr2];
    end
    
    mdl_cdf=zeros(51,numel(bin_edges));
    for z=1:51
        
        if layer==1
           idx=model_samples(z).smpl<0.1 | model_samples(z).smpl>0.9;
           model_samples(z).smpl(idx)=[];
        end
        
        f_sample=abs(model_samples(z).smpl-0.5)+0.5;
        Ni = sparse(histc(f_sample,bin_edges));
        li=sum(Ni);
        mdl_cdf(z,:)=cumsum(Ni)/li;
    end
    
    models_cdf(layer).models=mdl_cdf;
    model_05(layer).model=model_samples(1).smpl;
    all_models(layer).model=model_samples;
end
