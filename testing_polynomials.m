%% Test Section

%% 1.Polynomials

%% 1.1 degree = 15

clear all; close all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) x(1)^10+3*x(1)^6*x(2)^15+x(2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(40,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

fig = figure();

%subplot(4,3,1);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 1.2 degree = 25

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) x(1)^25+3*x(1)^6*x(2)^15+x(2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(40,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

%subplot(4,3,2);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 1.3 degree = 25 n = 20

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) x(1)^25+3*x(1)^6*x(2)^15+x(2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(40,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

%subplot(4,3,3);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 1.4 degree = 3

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) x(1)^3 + x(2)^2 + x(1);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(20,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value


%subplot(4,3,4);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 2.nonpolynomial functions

%% 2.1 sin,cos

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) sin(x(1))*cos(x(2));

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(40,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value

%subplot(4,3,5);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 2.1 exponential

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) exp(x(1))*exp(1-x(2));

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(40,g)
rel_error_1=zeros(n,1);value=zeros(n,1);step_error=zeros(n-1,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs((I-value(k))/I);
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value


%subplot(4,3,6);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on




