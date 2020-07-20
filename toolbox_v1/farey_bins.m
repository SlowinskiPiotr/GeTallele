function [bin_edges, bin_centres]=farey_bins(n)    

bin_centres=farey_sequence(n);

for_edges = bin_centres(:)';

binwidth = diff(for_edges);
for_edges = [0, for_edges(1:end-1)+binwidth/2, 1];
for_edges = full(real(for_edges));

% Shift bins so the interval is ( ] instead of [ ) for

bin_edges = for_edges + eps(for_edges);
bin_edges(1) = 0;
bin_edges(end) = 1;