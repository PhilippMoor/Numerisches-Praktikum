clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining smooth function f
f=@(x) exp(x(1));

% integrand over unit square
g_1=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1;x(2)]))*sqrt(det(m_tau'*m_tau))*x(1);


n=10;
I=Gauss_Quadrature(20,g_1)+Gauss_Quadrature(20,g_2)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_1(k)=abs(I-value(k));
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

fig = figure();

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on


%%

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) exp(x(1))*log(norm([x(1);x(2)]-A(1:2,1)));

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=10;
I=Gauss_Quadrature(20,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs(I-value(k));
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

fig = figure();

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on



