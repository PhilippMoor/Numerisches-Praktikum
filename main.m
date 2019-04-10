%% Main Programm
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

n=10;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(a) x^.2+3*x*y+y;

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

Gauss_Quadrature(n,g)

