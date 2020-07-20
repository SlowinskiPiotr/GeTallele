function [models_cdf,model_samples]=generate_VAF(totalreads,bin_edges,max_VAF)

for k=100:-1:50
    nb_reads=numel(totalreads);
    
    p=k/100;
    p=p*ones(1,5000);
    new_count_idx=randi(nb_reads,1,5000);
    new_count=totalreads(new_count_idx);
    y = binornd(new_count(:),p(:));
    surr_reads1=y(:)./new_count(:);
    
    p=1-k/100;
    p=p*ones(1,5000);
    new_count_idx=randi(nb_reads,1,5000);
    new_count=totalreads(new_count_idx);
    y = binornd(new_count(:),p(:));
    surr_reads2=y(:)./new_count(:);
    
    sample=[surr_reads1; surr_reads2];
    
    if max_VAF<1
        idx=sample>max_VAF | sample<(1-max_VAF);
        sample(idx)=[];
    end     

    model_samples(k-49).smpl=sample;
end

models_cdf=zeros(51,numel(bin_edges));

for z=51:-1:1
    sample=abs(model_samples(z).smpl-0.5)+0.5;
    Ni = sparse(histc(sample,bin_edges));
    li=sum(Ni);
    models_cdf(z,:)=cumsum(Ni)/li;
end

