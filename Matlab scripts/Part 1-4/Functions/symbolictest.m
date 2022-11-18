clear all
close all
clc

syms A [1 3]
B = 1+A(1)+A(2);
x = [1:3];

B_ans = subs(B,A,x)

x = 10*[1:3];
subs(B,A,x)

%%
x = 5*[1:3];
r = matlabFunction(subs(B,A,x), 'vars', 'x')
r(1)


%%
clear all
close all
clc

N = 10;
syms p
x =sym('x')

Ad = [1 0.1; -0.1-p 1.2]
% x = 10*[1:N];

A_quasi = matlabFunction(subs(Ad,p,x), 'vars', x)
x_in = 10*[1:N];
A_quasi(0)

Phi = 
