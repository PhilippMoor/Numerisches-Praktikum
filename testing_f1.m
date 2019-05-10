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

subplot(4,3,1);

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

subplot(4,3,2);

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

subplot(4,3,3);

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
    rel_error_1(k)=abs(I-value(k));
    if k>1
        step_error(k-1)=abs(value(k)-value(k-1));
    end
end
rel_error_1
step_error
value


subplot(4,3,4);

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

subplot(4,3,5);

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


subplot(4,3,6);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on



%% 3.Singularities

%% 3.1 singularity at Vertex of Triangle A type 1;

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
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

% subplot(4,3,7);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on


%% 3.1 singularity at Vertex of Triangle A type 2;

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=2;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
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

subplot(4,3,8);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 3.1 singularity at Vertex of Triangle A type 3;

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=3;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
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

subplot(4,3,9);

semilogy(1:n,rel_error_1,'blue*',1:n-1, step_error,'redo')
legend('relative error','step error')
xlabel('n');
grid on

%% 3.1 singularity at Vertex of Triangle A all types;

clear all; clc;

A=[1;2;0];B=[-1;0;3];C=[4;1;1];
type=1;

% Computing transformation matrix on unit triangle
m_tau=make_m_tau_mat(A,B,C,type);

% defining integrand over triangle
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);

% transforming to integrand over unit square
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);

n=20;
I=Gauss_Quadrature(20,g)
rel_error_1=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_1(k)=abs(I-value(k));
end
type=2;

m_tau=make_m_tau_mat(A,B,C,type);
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
n=20;
I=Gauss_Quadrature(20,g)
rel_error_2=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_2(k)=abs(I-value(k));
end

type=3;
m_tau=make_m_tau_mat(A,B,C,type);
f=@(x) 1/sqrt((x(1)-A(1))^2+(x(2)-A(2))^2);
g=@(x) f(chi(A,B,C,m_tau,type,rho(x)))*sqrt(det(m_tau'*m_tau))*x(1);
n=20;
I=Gauss_Quadrature(20,g)
rel_error_3=zeros(n,1);value=zeros(n,1);
for k=1:n
    value(k)=Gauss_Quadrature(k,g);
    rel_error_3(k)=abs(I-value(k));
end

hold off

% subplot(4,3,[10,11,12]);


semilogy(1:n,rel_error_1,'*',1:n,rel_error_2,'*',1:n,rel_error_3,'*')

legend('type 1','type 2','type 3')
xlabel('n');
grid on

h = figure(1);

print(h,'filename','-dpdf','-r0')


