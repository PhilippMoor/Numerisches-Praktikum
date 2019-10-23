%% Main Programm
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];

max=20;
error=zeros(3,max);


% defining integrand over triangle
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);

for type=1:3
    
    for n=1:10
        m_tau=make_m_tau_mat(A,B,C,type);
        g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
        Gauss_Quadrature(n,g)
    end
end

