function data=import_from_tsv(filename,purity)

from_file=importdata(filename);
nb_rows=size(from_file.textdata,1)-1;

%label
sampleID=from_file.textdata{2,5}(14:29);
our_label=double(string(from_file.textdata{2,5}(1:3)));

%layer order
lyr=NaN(1,4);

for i=1:4
    lr_str=from_file.textdata{1+i,5}(5:6);
    if strcmp(lr_str,'Ne')
        lyr(1)=i;
    elseif strcmp(lr_str,'Nt')
        lyr(2)=i;
    elseif strcmp(lr_str,'Te')
        lyr(3)=i;
    elseif strcmp(lr_str,'Tt')
        lyr(4)=i;
    else
        disp('Warning')
    end
end

for lr=1:4
    for i=lyr(lr):4:nb_rows
        lr_str=from_file.textdata{1+i,5}(5:6);
        if ~strcmp(lr_str,'Ne') && lr==1
            disp(filename)
            disp('wrong Nex layer - WHAT??? - in row:')
            disp(i)
            data=[];
            disp('correct the file and read data again')
            return
        end
        
        if ~strcmp(lr_str,'Nt')  && lr==2
            disp(filename)
            disp('wrong Ntr layer - WHAT??? - in row:')
            disp(i)
            data=[];
            disp('correct the file and read data again')
            return
        end
        
        if ~strcmp(lr_str,'Te')  && lr==3
            disp(filename)
            disp('wrong Tex layer - WHAT??? - in row:')
            disp(i)
            data=[];
            disp('correct the file and read data again')
            return
        end
        
        if ~strcmp(lr_str,'Tt')  && lr==4
            disp(filename)
            disp('wrong Ttr layer - WHAT??? - in row:')
            disp(i)
            data=[];
            disp('correct the file and read data again')
            return
        end
    end
end


%chromosome
cells4=from_file.textdata(2:4:end,1);
chromosomes=cellfun(@double,cellfun(@string,cells4,'UniformOutput',0));
chromosomes(isnan(chromosomes))=cellfun(@double,cellfun(@char,cells4(isnan(chromosomes)),'UniformOutput',0));
% 88 is X; 89 is Y; 77 is M

if mod(nb_rows,4)~=0
    pos=cellfun(@(x) str2double(x),from_file.textdata(2:end,2));
    for i=1:4:numel(pos)
        if 3*pos(i)-pos(i+1)-pos(i+2)-pos(i+3)~=0
            disp(filename)
            disp('too many or not enough layers - WHAT??? - near positions:')
            disp(pos(i:i+3))
            data=[];
            disp('correct the file and read data again')
            return
        end
    end
end

%position
positions=cellfun(@double,cellfun(@string,from_file.textdata(2:4:end,2),'UniformOutput',0));

%SNVcount 5
SNVcount=reshape(from_file.data(:,5),4,nb_rows/4)';
SNVcount=SNVcount(:,lyr);

%Refcount 6
Refcount=reshape(from_file.data(:,6),4,nb_rows/4)';
Refcount=Refcount(:,lyr);

%ratio 9
ratio=reshape(from_file.data(:,9),4,nb_rows/4)';
ratio=ratio(:,lyr);

data_matrix=[chromosomes, positions, ratio, SNVcount, Refcount];

idx_0109=(data_matrix(:,3)>=0.1 & data_matrix(:,3)<=0.9);
data_matrix_0109=data_matrix(idx_0109,:);

%purity
fun = @(s)any(strncmpi(s,sampleID,12));
purity_idx=cellfun(fun,purity.sample);
disp(['Number of purity estimates: ' num2str(sum(purity_idx))])
%out
data.purity=purity.values(purity_idx,:);
data.purity=data.purity(1,:);
data.matrix_raw=data_matrix;
data.matrix=data_matrix_0109;
data.sampleID=sampleID;
data.our_label=our_label;