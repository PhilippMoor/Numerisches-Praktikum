function m_tau = make_m_tau_mat(A,B,C,type)

% constructing m_tau (transforming matrix)
m_tau = zeros(3,2);

if type == 1        % -> A to the origin
    m_tau(:,1) = B-A;
    m_tau(:,2) = C-B;
elseif type == 2    % -> B
    m_tau(:,1) = C-B;
    m_tau(:,2) = A-C;
elseif type == 3    % -> C
    m_tau(:,1) = A-C;
    m_tau(:,2) = B-A;
else
    error('type must be 1,2 or 3');
end

end

