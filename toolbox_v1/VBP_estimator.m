function [VBP_est_keep,mix_complexity_keep,VBP_dist,all_prp]=VBP_estimator(chr_results,all_vprs,populations)
% This function estimates purity of a sample in the dataset. It uses
% matched normal and tumor exome signals (Nex and Tex). As returned by functions
% find_segments_and_fit_vpr run on Tex signal and later fit_vpr_in_segments
% run on Nex signal.
%
% Input: 
%   chr_results - structure with segmnts and estimatd vpr values for 
%                 22 chromosomes. Each of the 22 elements of the structure
%                 should have fields Tex and Nex; as outputed by functions:
%                 find_segments_and_fit_vpr run on Tex signal,
%                 fit_vpr_in_segments run on Nex signal
%                 run on that chromosomes
%   all_vprs - structure with 3D matrices containing all the possible vprs
%              values; output of gen_all_possible_vprs function.
%   populations - structure with all the possible population proportions;
%                 output of gen_all_possible_vprs function.
%
% Output:
%   VBP_est_keep - vpr based purity estimate
%   mix_complexity_keep - number of populations and events that ideentifies 
%                         the lowest complexity mixture (the 'simplest' mixture - 
%                         with the smallest number of populations and lowest
%                         maximal number of events).
%   VBP_dist - structure (number of populations x number of events) with all the possible purity values (1-proportion of 
%              normal population. The field .pr has all the admissible mixtures
%              for for a given number of populations and events
%   all_prp -  matrix with all the proportions (for any number of
%              populations and events) its columns are:  
%                       1. number of populations considered
%                       2. number of events considered
%                       3. proportion of the normal population
%                       4 to 5 (or 6). proportions of the tumor populations
%                       (depends on the number of populations provided for
%                       consideration in inputs all_vprs,populations)
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.      
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------


% collate all vprs from all the segments in all the chromosomes

all_vpr_Nex=NaN(1,22000); % preallocate memory assuming 1000 segments per chromosomes
all_vpr_Tex=NaN(1,22000); % preallocate memory assuming 1000 segments per chromosomes
all_Tex_p_05=NaN(1,22000); % preallocate memory assuming 1000 segments per chromosomes
all_seg_dp=NaN(1,22000); % preallocate memory assuming 1000 segments per chromosomes

% loop over all chromosomes
k=1;
for chr_idx=1:22
    data_and_stats_Tex=chr_results(chr_idx).Tex; % based on Tex
    data_and_stats_Nex=chr_results(chr_idx).Nex; % fitted in Nex
    
    vpr_Nex=data_and_stats_Nex.fitted_vprs; %copy estimated vprs values to new variable  
    vpr_Tex=data_and_stats_Tex.fitted_vprs; %copy estimated vprs values to new variable  
    p_05=data_and_stats_Tex.pv_vpr_neq_05; %copy p-values vallues to new variable  
    seg_dp=diff(data_and_stats_Tex.segment_edges); %compute lengths of the segments
    
    nb_of_segs=numel(vpr_Nex); %compute number of segments in the chromosome
    
    for i_segs=1:nb_of_segs 
        % in this loop we copy valuse from all segments in chromosome to
        % variabels that will have all the values from all the segments in
        % all chromosomes
        all_vpr_Nex(k)=vpr_Nex(i_segs);
        all_vpr_Tex(k)=vpr_Tex(i_segs);
        all_Tex_p_05(k)=p_05(i_segs);    
        all_seg_dp(k)=seg_dp(i_segs);
        k=k+1;
    end
end

% remove NaN values from the variables
all_vpr_Nex(isnan(all_vpr_Nex))=[];
all_vpr_Tex(isnan(all_vpr_Tex))=[];
all_Tex_p_05(isnan(all_Tex_p_05))=[];
all_seg_dp(isnan(all_seg_dp))=[];

% ------------- QUALITY CRITERIA FOR INCLUDING THE VPR VALUE IN THE PURITY ESTIMATION
% copy the vpr estimates from Tex signal to a new variable
vprs_for_VBP=all_vpr_Tex; 
% remove the vpr values for which vpr in the same segment in Nex signal was>0.58
vprs_for_VBP(all_vpr_Nex>0.58)=NaN; 
% remove the vpr values for which vpr in the same segment in Nex signal
% was no diffrent from synthethic VAF sample with vpr=0.5
vprs_for_VBP(all_Tex_p_05>1e-6)=NaN; 
% remove the vpr values estimated on less then 50 data points
vprs_for_VBP(all_seg_dp<50)=NaN;

% remove NaN values from the variables
vprs_for_VBP(isnan(vprs_for_VBP))=[];
% remove the vpr values <0.58
vprs_for_VBP(vprs_for_VBP<0.58)=[];

% ------------- FINDING ALL ADMISSIBLE MIXTURES
% check what is maximal number of population and events for which we will
% be serching for admissible mixtures
n_pop=size(all_vprs,1);
n_ev=size(all_vprs,2);

if ~isempty(vprs_for_VBP) %only start it if there are some good enough vpr values 
    % loop to find admissible mixtures
    for i_pop=n_pop:-1:2 %loop over populations
        for i_ev=n_ev:-1:1 %loop over events
            nb_prop=size(all_vprs(i_pop,i_ev).vprs,3);
            for prop=1:nb_prop %loop over considered mixtures
                vprs_p_e=all_vprs(i_pop,i_ev).vprs(:,:,prop);
                vprs_p_e=vprs_p_e(:); %all the vprs of a given mixture as a vector
                vprs_p_e(vprs_p_e<0.5)=[]; %folded remove all vpr values <0.5
                vprs_p_e=unique(vprs_p_e(:)); %unique vprs we only want uniqe values
                
                % find differences between all the vpr values of a synthetic
                % mixture and the observed/ estimated vpr values
                [~,dist]=knnsearch(vprs_p_e,vprs_for_VBP'); 
                
                % change status of the indicator variable that shows which
                % mixtures are admissible. if any of the observed vprs 
                % cannot be found among the synthethic mixture vprs. The
                % mixture is not c.
                
                if any(dist>0.009) || isempty(dist) %search radius is set 0.009
                    in_out(i_pop,i_ev).prp(prop)=NaN;   
                else
                    in_out(i_pop,i_ev).prp(prop)=1;
                end
            end
        end
    end
    
    % here we compute the total number of possible admissible mixtures to
    % pre-allocate output variables
    
    prp_times_ev_size=0;
    for i_pop=2:n_pop %loop over populations
        nb_prop=size(all_vprs(i_pop,2).vprs,3);
        prp_times_ev_size=prp_times_ev_size+n_ev*nb_prop;
    end
    
    % pre-allocating output variables
    prp_vec=zeros(prp_times_ev_size,n_pop);
    ev_vec=zeros(prp_times_ev_size,1);
    pop_vec=zeros(prp_times_ev_size,1);
    
    % here we read out of the indicator variable and copy, population
    % numbers, maximal event number and proportions to vectors and matrices
    k=1;
    for i_pop=2:n_pop %loop over populations
        for i_ev=1:n_ev %loop over events
            nb_prop=size(all_vprs(i_pop,i_ev).vprs,3);
            for prop=1:nb_prop %loop over considered mixtures
                if isnan(in_out(i_pop,i_ev).prp(prop))
                    prp_vec(k,:)=NaN;
                    ev_vec(k)=NaN;
                    pop_vec(k)=NaN;
                else
                    prp_vec(k,1:i_pop)=populations(1,i_pop).p(prop,:);
                    ev_vec(k)=i_ev;
                    pop_vec(k)=i_pop;
                end
                k=k+1;
            end
        end
    end
    
   
    all_prp=[pop_vec ev_vec prp_vec];  % big matrix with proportions of all admissible mixtures
    
    % here we convert proportion into purity=(1 - proportion of the normal population)
    % for all the admissible mixtures
    for p_vec=n_pop:-1:2
        for e_vec=n_ev:-1:1
            if sum(ev_vec==e_vec & pop_vec==p_vec)>0
                idx=(ev_vec==e_vec & pop_vec==p_vec); % we only consider admissible mixtures for a given number of proportions and events
                pr=round(1-prp_vec(idx,1),2);
                VBP_dist(p_vec,e_vec).pr=pr;
            else
                VBP_dist(p_vec,e_vec).pr=[];
            end
        end
    end
    
else
    % if there are no admissible mixtures return empty variables 
    VBP_dist=[];
    all_prp=[];
end

% ------------- SELECTION OF THE MIXTURE WITH THE LOWEST COMPLEXITY
% assignin output variables, after finding mixture with the lowest complexity 
if ~isempty(VBP_dist)
    to_keep=0;
    for p_i=2:n_pop % we 1st go over number of populations 
        for e_i=1:n_ev  % next we go over number of events
            if ~isempty(VBP_dist(p_i,e_i).pr)
                mix_complexity=[p_i,e_i]; % indicators of complexity (number of populations and events)
                VBP_est=min(VBP_dist(p_i,e_i).pr);
                disp('mixture found!')
                disp(['populations: ' num2str(p_i) '; events: ' num2str(e_i)]);
                disp(['Purity: ' num2str(VBP_est)]);
                
                % we mark the 1st found solution (it will have the smallest
                % number of populations, but might have high number of
                % events) as the one with the lowest complexity
                if to_keep==0 
                    disp('this mixture has the lowest complexity and this value will be returned');
                    mix_complexity_keep=mix_complexity;
                    VBP_est_keep=VBP_est;
                    to_keep=1; %this if statement will be entered only once
                end
            else
                disp(['no fit - populations: ' num2str(p_i) '; events: ' num2str(e_i)])
                if to_keep==0 % we only want NaN if there wasn't any admissible mixture for all populations and events
                    mix_complexity_keep=NaN;
                    VBP_est_keep=NaN;
                end

            end
        end
    end
else
    % if there are no vpr values left after applying quality criteria return empty variables 
    disp('no vprs')
    mix_complexity_keep=NaN;
    VBP_est_keep=NaN;
end