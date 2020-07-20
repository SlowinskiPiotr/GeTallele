[~,~,purity]=xlsread('./ncomms9971-s2.xlsx','Supp Data 1');

idx_B=1;

for i=5:size(purity,1)
    if strcmp(purity{i,2},'BRCA')
        purity_BRCA.sample{idx_B}=purity{i,1};
        for k=1:5
            if purity{i,k+2}=='NaN'
                purity_BRCA.values(idx_B,k)=NaN;
            else
                purity_BRCA.values(idx_B,k)=purity{i,k+2};
            end
        end
        idx_B=idx_B+1, 
    end
end


