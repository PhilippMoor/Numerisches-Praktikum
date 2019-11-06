%% testing one type
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining smooth function f
f=@(x) exp(x(1)^2+x(2)^2);

% integrand over unit square
g_1=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1;x(2)]))*sqrt(det(m_tau'*m_tau))*x(1);


n=20;
I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(k,g_2);
    rel_error_1(k)=abs((I-value(k))/I);
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

%% testing all types
clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
P=B;n=20;

% defining r
r=@(x) exp(x(1)^2+x(2)^2);

% defining f
f=@(x) r(x)*log(sqrt((x(1)-P(1))^2+(x(2)-P(2))^2));

type=1;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)
rel_error_1=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(20,g_2);
    rel_error_1(k)=abs((I-value(k))/I);
end

type=2;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)
rel_error_2=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(20,g_2);
    rel_error_2(k)=abs((I-value(k))/I);
end

type=3;
m_tau=make_m_tau_mat(A,B,C,type);
g_1=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(x(1))*sqrt(det(m_tau'*m_tau))*x(1);
g_2=@(x) r(chi(A,B,C,m_tau,type,rho(x)))*log(norm(m_tau*[1,x(2)]'))*sqrt(det(m_tau'*m_tau))*x(1);

I=Gauss_Quadrature(40,g_1)+Gauss_Quadrature(40,g_2)
rel_error_3=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g_1)+Gauss_Quadrature(20,g_2);
    rel_error_3(k)=abs((I-value(k))/I);
end

% subplot(4,3,[10,11,12]);

h = figure();
semilogy(1:n,rel_error_1,'*',1:n,rel_error_2,'*',1:n,rel_error_3,'*')

legend('type A','type B','type C')
xlabel('n');
grid on
%saveas(h,'g_(x,y)','png');
