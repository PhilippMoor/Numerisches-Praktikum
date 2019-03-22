clear all;close all;clc;

% dimension n = #number of nodes in one coordinate
n = 2;
% f = function to be integrated over the ->square<-

% load structure with nodes and weights
struct = load("nodes.mat");

x = struct.nodes(n+1,1:n+1);
y = struct.nodes(n,1:n);
weights_x = struct.weights(n+1,1:n+1);
weights_y = struct.weights(n,1:n);

% definition of rho
rho=@(a,b)[a,a*b];
