clear all;
close all;


m = 2:1:7; % bcakoff stage

syms M W0 I Xi p

Wi = W0*2^I;
Xi = (-1*M/2)*floor(Wi/M)^2 + (Wi-M/2)*floor(Wi/M);

for k = 1:max(m)
    k
    Wi = W0*2^k;
    X(k) = (-1*M/2)*floor(Wi/M)^2 + (Wi-M/2)*floor(Wi/M)
end

% for m = 2:7
%     m
%     symsum((p/2)^I, I, 1, m-1)
%     
% end