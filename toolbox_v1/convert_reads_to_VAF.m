function data_VAF=convert_reads_to_VAF(data_readcounts)
% This function computes VAF values from readcounts. It keeps the
% position (chromosome and base pair) for each VAF value
%
% Call function: data_VAF=convert_reads_to_VAF(data_readcounts)
%
%      Input: 
%            data_readcounts - matrix:
%                              1st column - chromosome
%                              2nd column - base pair on the chromosome
%                              3rd column - variant readcount
%                              4th column - reference readcount 
%      Output:
%            data_VAF - matrix:
%                              1st column - chromosome
%                              2nd column - base pair on the chromosome
%                              3rd column - VAF value
%        
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

data_VAF=[data_readcounts(:,[1 2]) data_readcounts(:,3)./(data_readcounts(:,3)+data_readcounts(:,4))];
