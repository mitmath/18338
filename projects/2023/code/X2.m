function [x] = X2(t,s)
%X1 Summary of this function goes here
%   Detailed explanation goes here
x=zeros([size(t),3]);
x(:,:,1)=zeros(size(t));
x(:,:,2)=s;
x(:,:,3)=2*t.*s;
end

