function [all_vprs,populations,alleles]=gen_all_possible_vprs(p,a,ad_mix_res)
% This function generates all possible vprs values that can be observed
% in a segment for a given number of populations (mixed in all possible
% ratios) and some maximal number of genetic events occuring (such as 
% deletions, duplications etc. all the way to pentaplications.
%
% Input: 
%   p - maximal number of populations (for memory reasons max(p)=4, but it could be changed in future)
%   a - maximal number of alleles' 'multiplicity' (for memory reasons max(a)=5, but it could be changed)
%       ----- changing the maximal allowed 'p' and 'a' will require changing lines 126-152 ----------------    
%   ad_mix_res - 'resolution' of considerd mixtures, e.g. 0.01 means that all
%                proportions of the populations will be considered from 0.01:0.01:0.99 will be considered
%
% Output:
%   all_vprs - structure with matrices containing all the possible vprs values, e.g.
%                all_vprs(4,3) contains vprs values for mixtures of 4
%                populations with up to 3 alleles, and all_vprs(3,2)
%                contains vprs values for mixtures of 3 populations with
%                up to 2 alleles. All the elements of the matrices are vpr values.
%                Matrices have 3 dimensions: ev_a x ev_b x mix, with
%                 ev_a - combinations of events allele a     
%                 ev_b - combinations of events allele b    
%                 mix - mixtures (different ratios of populations in the sample).
%
%   populations - structure with all the possible population proportions,
%               e.g. populations(2) contains all the considered proportions
%               of 2 populations
%   alleles - structure with all the alleles combinations for different
%           numbers of populations, e.g. alleles(4,3) contains 4 elements combinations
%           (4 populations) of 3 element set (3 alleles)
%
%   structures for 1 popuation are empty
%
% Function uses the combinator function by Matt Fig 
% (just Matlab without the mex option) which is include in the \external\ternplot folder.
% 
% Matt Fig (2020). COMBINATOR -combinations AND permutations 
% (https://www.mathworks.com/matlabcentral/fileexchange/24325-combinator-combinations-and-permutations), 
% MATLAB Central File Exchange. Retrieved July 24, 2020. 
%
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.      
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

if a>5
    warning('a - alleles multiplicity is to high')
    return
end

if p>4
    warning('to many populations')
    return
end

% make all possible combinations of allels in segments in the tumor population
for i_al=a:-1:1
    for i_pop=p:-1:2
        alleles(i_pop,i_al).a=combinator(i_al+1,i_pop-1,'c','r')-1; % COMBINATIONS WITH REPETITION
        alleles(i_pop,i_al).b=combinator(i_al+1,i_pop-1,'c','r')-1; % COMBINATIONS WITH REPETITION
        % matrices with numbers from 0 to 5
    end
end

% make all possible permutations of population ratios (for a given resolution)
for i=p:-1:2
    p_pre=(combinator(round(1/ad_mix_res+1),i-1,'p','r')-1)*ad_mix_res; % PERMUTATIONS WITH REPETITION
    % e.g for i=1 (single population) and ad_mix_res=0.01 p_pre is a vector
    % with numbers 0:0.01:1 (it starts as numbers 1:101 but we subtract 1
    % and multiply by 0.01)
    
    p_pre=round(p_pre,2); % making sure there only two decimal places
    p_pre=p_pre(sum(p_pre,2)<1,:); 
    % include only this sets of ratios that sum to <1
    % so we also exclude mixtures with one of the proportions 1
    p_all=[p_pre 1-sum(p_pre,2)]; 
    % add the last ratio (and now they sum to 1) 
    p_all=round(p_all,2); 
    % again making sure there only two decimal places

    i0=(p_all==0);  % find all the 0 ratios
    populations(i).p=p_all(sum(i0,2)==0,:); 
    % exclude any mixtur with a 0 ratio (it is included as a one of the 
    % mixtures with 1 less population)
end

% make matrices with all vprs values with for all mixtures and all events
% combinations
for i_al=a:-1:1 % loop over events
    for i_pop=p:-1:2 % loop over number of populations
        disp(['now computing all vprs for: ' num2str(i_al) ' events, ' num2str(i_pop) ' populations'])
        sz_alA=size(alleles(i_pop,i_al).a,1); % finding sizes to run loops and preallocate memory
        sz_alB=size(alleles(i_pop,i_al).b,1); % finding sizes to run loops and preallocate memory
        sz_pop=size(populations(i_pop).p,1); % finding sizes to run loops and preallocate memory
        all_vprs(i_pop,i_al).vprs=NaN(sz_alA,sz_alB,sz_pop); % make an empty matrix that will contain all the vpr values
        % the 3D matrix has dimennsion sz_alA x sz_alB x sz_pop, 
        % sz_alA - number of possible combinations of events of allele a     
        % sz_alB - number of possible combinations of events of allele b    
        % sz_pop - number of all possible mixtures.
        
        % this loop populates the matrix with all possible vpr values
        for iA=sz_alA:-1:1
            for iB=sz_alB:-1:1
                for i_pop_prop=1:sz_pop
                    % take combinations of events and multiply them by
                    % proportions in the mixture
                    % we add 1 here as the 1st elements to include 1:1 ratio 
                    % of the normal population
                    A=sum([1 alleles(i_pop,i_al).a(iA,:)].*populations(i_pop).p(i_pop_prop,:));
                    B=sum([1 alleles(i_pop,i_al).b(iB,:)].*populations(i_pop).p(i_pop_prop,:));
                    % compute probability of seeing one of the alleles 
                    all_vprs(i_pop,i_al).vprs(iA,iB,i_pop_prop)=B/(A+B);
                    % (we use all the poossible combinations of events so the symmetric counter parts are included)
                end
            end
        end
        
% finally, we filter out vprs values that would be produced if one of the population
% had a segment completly removed. In principle entries to the vprs
% matrix allow deletion of both alleles in a given segment in one or more 
% populations. But we consider such an event highly unlikely. 

% To remove this filtering comment this part of the code        
        if i_pop==2 % filter if 2 populations
            all_vprs(i_pop,i_al).vprs(1,1,:)=0;
        elseif i_pop==3 % filter if 3 populations
            if i_al==1
                all_vprs(i_pop,i_al).vprs(1:2,1:2,:)=0;
            elseif i_al==2
                all_vprs(i_pop,i_al).vprs(1:3,1:3,:)=0;
            elseif i_al==3
                all_vprs(i_pop,i_al).vprs(1:4,1:4,:)=0;
            elseif i_al==4
                all_vprs(i_pop,i_al).vprs(1:5,1:5,:)=0;
            elseif i_al==5
                all_vprs(i_pop,i_al).vprs(1:6,1:6,:)=0;
            end
        elseif i_pop==4 % filter if 4 populations
            if i_al==1
                all_vprs(i_pop,i_al).vprs(1:3,1:3,:)=0;
            elseif i_al==2
                all_vprs(i_pop,i_al).vprs(1:6,1:6,:)=0;
            elseif i_al==3
                all_vprs(i_pop,i_al).vprs(1:10,1:10,:)=0;
            elseif i_al==4
                all_vprs(i_pop,i_al).vprs(1:15,1:15,:)=0;
            elseif i_al==5
                all_vprs(i_pop,i_al).vprs(1:21,1:21,:)=0;
            end
        end        
     end
end
