function [m_tau,determinant] = make_m_tau_mat(A,B,C,fixed_vertex)
% convert 3D input into 2D

A = A(1:2)';
B = B(1:2)';
C = C(1:2)';

fixed_vertex = fixed_vertex';

% constructing m_tau (transforming matrix)

m_tau = zeros(2,2);
m_tau(:,1) = B-A;
m_tau(:,2) = C-B;

% determinant (D_chi_tau)



end

