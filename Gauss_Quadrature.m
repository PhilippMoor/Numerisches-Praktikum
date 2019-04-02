clear all;close all;clc;

% dimension n = #number of nodes in one coordinate
n = 10;
% f = function to be integrated over the ->square<-

% load structure with nodes and weights
struct = load("nodes.mat");

x = struct.nodes(n+1,1:n+1);
y = struct.nodes(n,1:n);
weights_x = struct.weights(n+1,1:n+1);
weights_y = struct.weights(n,1:n);

% definition of rho
rho_1=@(a,b)a;
rho_2=@(a,b)a*b;

%definition of f
f=@(x,y) x.^2+y.^2;

% function after transformation
f_rho=@(a,b) f(rho_1(a,b),rho_2(a,b))*a;

% Computing the Quadrature
sum=0;
for i=1:n+1
    sum_y=0;
    for j=1:n
        sum_y=sum_y+(weights_x(i)*weights_y(j)*f_rho(x(i),y(j)));
    end
    sum=sum+sum_y;
end
disp(sum)

