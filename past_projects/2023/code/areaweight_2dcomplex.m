function [func] = areaweight_2dcomplex(t,s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c=1 + t.^2 + (t.^2.-s.^2).^2;
d=1+4.*t.^2+(t.^2-s.^2).^2;
f=1+4*t.^2;
func=2*s.^2.*(d.^-0.5).*abs((f.^-0.5) - s.^2.*(c.^-0.5));
end

