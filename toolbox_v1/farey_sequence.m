function f_ss=farey_sequence(n)
%function [f_s, ford_r]=farey_sequence(n)

a=0;
b=1;
c=1;
d=n;

i=1;
f_ss=NaN(1,ceil((3*n^2)/pi^2+10000));
f_s=0;

while f_s<1 
    p=a;
    q=b;
    f_s=p/q;
    f_ss(i)=f_s;
    %ford_r(i)=1/(2*q^2);
    i=i+1;
    
    k=floor((n+b)/d);
    a=c;
    b=d;
    c=(k*c-p);
    d=(k*d-q);
end

f_ss(isnan(f_ss))=[];