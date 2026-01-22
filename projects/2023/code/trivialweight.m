function [func] = trivialweight(t,s)
%TRIVIALWEIGHT Summary of this function goes here
%   Detailed explanation goes here
func=exp(-t.^2-s.^2);
end

