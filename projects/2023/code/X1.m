function [x] = X1(t,s)
%X1 Summary of this function goes here
%   Detailed explanation goes here
x=zeros([size(t),3]);
x(:,:,1)=ones(size(t));
x(:,:,2)=t;
x(:,:,3)=t.^2-s.^2;
end

