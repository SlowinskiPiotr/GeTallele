function convert_VAF_to_circos_data(data_mat,fname)

fid=fopen(fname,'w');

for chrm=1:22
    chr_idx=data_mat(:,1)==chrm;
    chr_data=data_mat(chr_idx,2:3);
    for data_pts=1:size(chr_data,1)
        pos=chr_data(data_pts,1);
        value=chr_data(data_pts,2);
        datastr=['hs' num2str(chrm) ' ' num2str(pos) ' ' num2str(pos) ' ' num2str(value) '\n'];
        fprintf(fid,datastr);
    end
end
fclose(fid);