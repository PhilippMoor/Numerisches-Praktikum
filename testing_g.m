clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
P=A;
n=20;

% defining f
f=@(x) log(sqrt((x(1)-P(1))^2+(x(2)-P(2))^2));

type=1;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
I=Gauss_Quadrature(20,g)
rel_error_1=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
end

type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
I=Gauss_Quadrature(20,g)
rel_error_2=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_2(k)=abs((I-value(k))/I);
end

type=3;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
I=Gauss_Quadrature(20,g)
rel_error_3=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_3(k)=abs((I-value(k))/I);
end


h = figure();
semilogy(1:n,rel_error_1,'*',1:n,rel_error_2,'*',1:n,rel_error_3,'*')
legend('type A','type B','type C')
xlabel('n');
grid on
