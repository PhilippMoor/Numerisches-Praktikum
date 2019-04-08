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

%% Test Section
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) x(1)^10+3*x(1)^6*x(2)^15+x(2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=10;
I=Gauss_Quadrature(20,g)
rel_error=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error(k)=abs(I-value(k));
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error
step_error
value

plot(1:n,rel_error,'blue*',1:n-1, step_error,'redo')

semilogy(1:n,value)
grid on


