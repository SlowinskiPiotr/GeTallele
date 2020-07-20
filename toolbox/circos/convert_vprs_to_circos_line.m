function convert_vprs_to_circos_line(data_mat,layer,loc_half,loc_width,top_or_bottom,fname)

fid=fopen(fname,'w');

for chrm=1:22
    chr_idx=data_mat(:,1)==chrm;
    chr_data=data_mat(chr_idx,[2:3 4:7 8:11 12 17]); %1:12
    for data_pts=1:size(chr_data,1)
        start_pos=chr_data(data_pts,1);
        end_pos=chr_data(data_pts,2);
        value=chr_data(data_pts,layer+2);
        
        if layer==1 || layer==2
            pv05l1=chr_data(data_pts,7)>1e-5; %not different from 0.5
            pv05l2=chr_data(data_pts,8)>1e-5; %not different from 0.5
            pvl1l2=chr_data(data_pts,11)>1e-5; %same
            
            ex_c1=(pv05l1~=pv05l2) & pvl1l2; % one differtn from 0.5, one not different from 0.5 but same
        elseif layer==3 || layer==4
            pv05l1=chr_data(data_pts,9)>1e-5;  %not different from 0.5
            pv05l2=chr_data(data_pts,10)>1e-5; %not different from 0.5
            pvl1l2=chr_data(data_pts,12)>1e-5; %same
            
            ex_c1=(pv05l1~=pv05l2) & pvl1l2; % one differtn from 0.5, one not different from 0.5 but same
        end
        
        if layer==1 || layer==3
            ex_c2=chr_data(data_pts,layer+6)>1e-5 & value>=0.58; %not different from 0.5 but have high vpr
        elseif layer==2 || layer==4
            ex_c2=chr_data(data_pts,layer+6)>1e-5 & value>=0.58; %not different from 0.5 but have high vpr
        end
        
        if ~(ex_c1 || ex_c2)
            if chr_data(data_pts,layer+6)>1e-5
                value=0.5;
            end
            if (layer==1 || layer==3) && value<0.58
                value=0.5;
            end
            if (layer==2 || layer==4) && value<0.58
                value=0.5;
            end
            
            if top_or_bottom==1
                value=loc_half+loc_width*(value-0.5)/0.5;
            else
                value=loc_half-loc_width*(value-0.5)/0.5;
            end
            
            datastr=['hs' num2str(chrm) ' ' num2str(start_pos) ' ' num2str(end_pos) ' ' num2str(value) ' ' ...
                'r0=' num2str(value,5) 'r-3p,' 'r1=' num2str(value,5) 'r+3p' '\n'];
            fprintf(fid,datastr);
        end
    end
end
fclose(fid);
