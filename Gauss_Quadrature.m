clear all;close all;clc;

% dimension n = #number of nodes in one coordinate
n = 2;
% f = function to be integrated over the ->square<-

% load matrices with nodes and weights
nodes = load("nodes.mat")

x = nodes.nodes(n,1:n+1)
y = x
weights_x = nodes.weights(n,1:n)
weights_y = weights_x
