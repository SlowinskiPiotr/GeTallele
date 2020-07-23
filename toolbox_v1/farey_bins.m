function [bin_edges, bin_centres]=farey_bins(n)    
% This function generates bin edges, centered at the elements of the Faray sequence of order n.

% Call function: f_seq=farey_sequence(n)
%
%      Input: 
%            n - order of the Farey sequence
%      Output:
%            bin_edges - vector of bin edges based on Farey sequence of order n, with number of elements equal approximatley to (3*n^2)/pi^2 
%            bin_centres - vector with Farey sequence of order n, with number of elements equal approximatley to (3*n^2)/pi^2 

% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

bin_centres=farey_sequence(n); %compute faree sequence of order n

binwidth = diff(bin_centres); %compute widths of the individual bins
for_edges = [0, bin_centres(1:end-1)+binwidth/2, 1]; %keep the limits of 0 and 1, add half of the width to get all the bin edges in between 0 and 1

% Shift bins so the interval is ( ] instead of [ )
bin_edges = for_edges + eps(for_edges);

% again make sure the the limits are 0 and 1
bin_edges(1) = 0; 
bin_edges(end) = 1;