function [integral,error] = Gauss_Quadrature(n,f)

% load structure with nodes and weights
struct = load("nodes.mat");

x = struct.nodes(n+1,1:n+1);
y = struct.nodes(n,1:n);
weights_x = struct.weights(n+1,1:n+1);
weights_y = struct.weights(n,1:n);



% Computing the quadrature
sum=0;
for i=1:n+1
    sum_y=0;
    for j=1:n
        nodes=[x(i),y(j)];
        sum_y=sum_y+(weights_x(i)*weights_y(j)*f(nodes));
    end
    sum=sum+sum_y;
end
integral=sum;


% Computing the error
error=1;


end