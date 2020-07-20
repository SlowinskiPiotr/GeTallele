[~,~,genecoordinates]=xlsread('./hg38_gene_coordinates.xlsx','genes');

nb_genes=size(genecoordinates,1);
gene_names={genecoordinates{2:end,2}};

chr_pos=NaN(nb_genes-1,3);

for i=2:nb_genes
    coord_stra=genecoordinates{i,3};
    idx_col=find(coord_stra==':');
    idx_bar=find(coord_stra=='-');
    chrom=string(coord_stra(4:(idx_col-1)));
    if strcmp(chrom,'M')
        chr_pos(i-1,1)=100;
    elseif strcmp(chrom,'X')
        chr_pos(i-1,1)=23;
    elseif strcmp(chrom,'Y')
        chr_pos(i-1,1)=24;
    else
        chr_pos(i-1,1)=double(chrom);
    end
    chr_pos(i-1,2)=double(string(coord_stra((idx_col+1):(idx_bar-1))));
    chr_pos(i-1,3)=double(string(coord_stra((idx_bar+1):end)));
end

save gene_coordinates gene_names chr_pos


