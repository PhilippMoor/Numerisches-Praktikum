%% testing one type
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;n=20;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining smooth function f
r=@(x) exp(x(1)^2+x(2)^2);

% integrand over unit square
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1;x(2)]))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_1(k)=abs(I-value(k))/I;
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

fig = figure();
semilogy(1:n,rel_error_1,'blue*')
legend('relative error','step error')
xlabel('n');
grid on

%% testing all types
clear all; close all; clc;

A=[1;2;0];B=[1;0;3];C=[4;1;1];
n=20;

% defining r
r=@(x) exp(x(1)^2+x(2)^2);

% Computing overkill
type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)

% Computing with different types
type=1;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_1=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_1(k)=abs(value(k)-I)/I;
end

type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_2=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_2(k)=abs(value(k)-I)/I;
end

type=3;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

rel_error_3=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_3(k)=abs(value(k)-I)/I;
end

fig = figure();
semilogy(1:n,rel_error_1,'*',1:n,rel_error_2,'*',1:n,rel_error_3,'*')
legend('type A','type B','type C')
xlabel('n');
grid on
