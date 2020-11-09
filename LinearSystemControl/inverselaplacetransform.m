clear all; close all; clc;
syms s;

A = [-1 -1 0 ; 0 -1 0 ; 0 0 -1];
% A = [-1 -1 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 -1]
ilaplace(inv(s*eye(size(A,1))-A))