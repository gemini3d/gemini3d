function [int,dx]=intrap(f,x)

%Integrates along a column using the trapezoidal rule
%
%function [int,dx]=intrap(f,x)
% bounds=size(f);
% dx=abs(diff(x,1,1));
% fmid=(f(1:bounds(1)-1,:)+f(2:bounds(1),:))/2;
% 
% int=cumsum(fmid.*dx,1);

bounds=size(f);
dx=abs(diff(x,1,1));
fmid=(f(1:bounds(1)-1,:)+f(2:bounds(1),:))/2;

int=cumsum(fmid.*dx,1);

end