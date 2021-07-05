clear all;
close all;

syms k W M

A = floor((W-k)/M) + 1;

Sn = symsum(A, k, M+1, W)