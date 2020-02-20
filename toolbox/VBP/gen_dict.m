function [dictionaries,populations,alleles]=gen_dict_old(p,a,ad_mix_res)
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
%
% p - maximal number of populations (for memory reasons max(p)=4 could be changed in future)
% a - maximal number of alleles' 'multiplicity' (for memory reasons max(a)=5 could to be changed in future)
%
% output:
% dictionaries - structure with all the possible dictionaries, e.g.
%                dictionaries(4,3) contains dictionary for mixtures of 4
%                populations with up to 3 alleles, and dictionaries(3,2)
%                contains dictionary for mixtures of 3 populations with
%                up to 2 alleles.
% populations - structure with all the possible population proportions,
%               e.g. populations(2) contains all the considered proportions
%               of 2 populations
% alleles - structure with all the alleles combinations for different
%           numbers of populations, e.g. alleles(4,3) contains 4 elements combinations
%           (4 populations) of 3 element set (3 alleles)
%
% structures for 1 popuation are empty
%
% If you have any questions please contact: p.m.slowinski@exeter.ac.uk


alleles=struct([]);
populations=struct([]);
dictionaries=struct([]);

if a>5
    warning('a - alleles multiplicity is to high')
    return
end

if p>4
    warning('to many populations')
    return
end

for i_al=1:a
    for i_pop=1:(p-1)
        alleles(i_pop+1,i_al).a=combinator(i_al+1,i_pop,'c','r')-1;
        alleles(i_pop+1,i_al).b=combinator(i_al+1,i_pop,'c','r')-1;
        %alleles(i_pop+1,i_al).b=combinator(2,i_pop,'c','r')-1;
    end
end


for i=2:p
    p_pre=(combinator(round(1/ad_mix_res+1),i-1,'p','r')-1)*ad_mix_res;
    %p_pre=(combinator(51,i-1,'p','r')-1)*2/100;
    %p_pre=(combinator(21,i-1,'p','r')-1)*5/100;
    p_pre=round(p_pre,2);
    p_pre=p_pre(sum(p_pre,2)<1,:);
    p_all=[p_pre 1-sum(p_pre,2)];
    p_all=round(p_all,2);
    p_all=p_all(2:end,:);
    i0=(p_all==0);
    populations(i).p=p_all(sum(i0,2)==0,:);
end

for i_al=1:a
    for i_pop=2:p
        tic,
        disp(['Dictionary: ' num2str(i_al) ' events, ' num2str(i_pop) ' populations'])
        sz_alA=size(alleles(i_pop,i_al).a,1);
        sz_alB=size(alleles(i_pop,i_al).b,1);
        sz_pop=size(populations(i_pop).p,1);
        dictionaries(i_pop,i_al).dict=NaN(sz_alA,sz_alB,sz_pop);
        
        for iA=1:sz_alA
            for iB=1:sz_alB
                for i_pop_prop=1:sz_pop
                    A=sum([1 alleles(i_pop,i_al).a(iA,:)].*populations(i_pop).p(i_pop_prop,:));
                    B=sum([1 alleles(i_pop,i_al).b(iB,:)].*populations(i_pop).p(i_pop_prop,:));
                    dictionaries(i_pop,i_al).dict(iA,iB,i_pop_prop)=B/(A+B);
                end
            end
        end
        
        if i_pop==2
            dictionaries(i_pop,i_al).dict(1,1,i_pop_prop)=0;
        elseif i_pop==3
            if i_al==1
                dictionaries(i_pop,i_al).dict(1:2,1:2,i_pop_prop)=0;
            elseif i_al==2
                dictionaries(i_pop,i_al).dict(1:3,1:3,i_pop_prop)=0;
            elseif i_al==3
                dictionaries(i_pop,i_al).dict(1:4,1:4,i_pop_prop)=0;
            elseif i_al==4
                dictionaries(i_pop,i_al).dict(1:5,1:5,i_pop_prop)=0;
            end
        elseif i_pop==4
            if i_al==1
                dictionaries(i_pop,i_al).dict(1:3,1:3,i_pop_prop)=0;
            elseif i_al==2
                dictionaries(i_pop,i_al).dict(1:6,1:6,i_pop_prop)=0;
            elseif i_al==3
                dictionaries(i_pop,i_al).dict(1:10,1:10,i_pop_prop)=0;
            elseif i_al==4
                dictionaries(i_pop,i_al).dict(1:15,1:15,i_pop_prop)=0;
            end
        end
        
    end
end
