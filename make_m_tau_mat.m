function m_tau = make_m_tau_mat(A,B,C,type)

% constructing m_tau (transforming matrix)
m_tau = zeros(3,2);

if type == 1
    m_tau(:,1) = B-A;
    m_tau(:,2) = C-B;
elseif type == 2
    m_tau(:,1) = C-B;
    m_tau(:,2) = A-C;
elseif type == 3
    m_tau(:,1) = A-C;
    m_tau(:,2) = B-A;
else
end

end

