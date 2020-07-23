function f_seq=farey_sequence(n)
% This function generates a Faray sequence of order n. The code is based on https://rosettacode.org/wiki/Farey_sequence
%
% Call function: f_seq=farey_sequence(n)
%
%      Input: 
%            n - order of the Farey sequence
%      Output:
%            f - vector with Farey sequence of order n, with number of elements equal approximatley to (3*n^2)/pi^2 
%
% The authors make no representations about the suitability of this software for any purpose.
% It is provided "as is" without express or implied warranty.
% -------------------------------------------------------------------------
%   P. Slowinski, p.m.slowinski@exeter.ac.uk, 2020
% -------------------------------------------------------------------------

% constants
a=0;
b=1;
c=1;
d=n;

i=1;

%pre-allocatio of variable
f_seq=NaN(1,ceil((3*n^2)/pi^2+10000)); %we add 10000 because we only have an estimate of the number of elements.

f_s=0; % the 1st  element is equal to 0

% loop runs until we get element of the Farey sequence equal to 1
while f_s<1 
    p=a;
    q=b;
    f_s=p/q; %compute the element
    f_seq(i)=f_s; %store the element
    i=i+1;
    
    k=floor((n+b)/d);
    a=c;
    b=d;
    c=(k*c-p);
    d=(k*d-q);
end

% remove unused elements
f_seq(isnan(f_seq))=[];