clear all; close all; clc;

A=[0;0;0];B=[1;0;0];C=[1;1;0];
P=B;n=20;

% defining f
f=@(x) log(sqrt((x(1)-P(1))^2+(x(2)-P(2))^2));

% Computing overkill
type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g)

% Computing with different types
type=1;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_1=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((value(k)-I)/I);
end

type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_2=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_2(k)=abs((value(k)-I)/I);
end

type=3;
m_tau=make_m_tau_mat(A,B,C,type);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_3=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_3(k)=abs((value(k)-I)/I);
end

fig = figure();
semilogy(1:n,rel_error_1,'o',1:n,rel_error_2,'*',1:n,rel_error_3,'*')
legend('type A','type B','type C')
xlabel('n');
grid on
