function [sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_LLL( ...
          sample_K,sample_X,sample_Y,sample_Z,rec,...
          b,params)

m_idx = find(b);
M_idx = length(m_idx);

if M_idx<=1
    return
end

%% transfer in
sample_X_idx = sample_X(:,m_idx);
sample_Y_idx = sample_Y(:,m_idx);
sample_Z_idx = sample_Z(:,m_idx);
sample_K_idx = sample_K(  m_idx);

%% re-link
A = params.D_prior_A+1.5*(params.N-1)*M_idx;
sample_B_idx = params.D_prior_B + 0.25*sum([diff(sample_X_idx),diff(sample_Y_idx),diff(sample_Z_idx)].^2 ./ params.dt,[1 2]);

for rep = 1 : poissrnd(2*(params.N-1)*M_idx*(M_idx-1))

    if rand<0.5
        n = 1 + randi(params.N-1);
        n_idx = n:params.N;
    else
        n =     randi(params.N-1);
        n_idx = 1:n;
    end
        
    m1 = randi(M_idx  );
    m2 = randi(M_idx-1);
    m2 = m2 + (m2>=m1);
    
    propos_X_idx = sample_X_idx;
    propos_Y_idx = sample_Y_idx;
    propos_Z_idx = sample_Z_idx;
    propos_X_idx(n_idx,[m1 m2]) = sample_X_idx(n_idx,[m2 m1]);
    propos_Y_idx(n_idx,[m1 m2]) = sample_Y_idx(n_idx,[m2 m1]);
    propos_Z_idx(n_idx,[m1 m2]) = sample_Z_idx(n_idx,[m2 m1]);
        
    if (propos_X_idx(sample_K_idx(m1),m1)<params.X_prior_min) || (propos_X_idx(sample_K_idx(m1),m1)>params.X_prior_max) ...
    || (propos_Y_idx(sample_K_idx(m1),m1)<params.Y_prior_min) || (propos_Y_idx(sample_K_idx(m1),m1)>params.Y_prior_max) ...
    || (propos_Z_idx(sample_K_idx(m1),m1)<params.Z_prior_min) || (propos_Z_idx(sample_K_idx(m1),m1)>params.Z_prior_max) ...
    || (propos_X_idx(sample_K_idx(m2),m2)<params.X_prior_min) || (propos_X_idx(sample_K_idx(m2),m2)>params.X_prior_max) ...
    || (propos_Y_idx(sample_K_idx(m2),m2)<params.Y_prior_min) || (propos_Y_idx(sample_K_idx(m2),m2)>params.Y_prior_max) ...
    || (propos_Z_idx(sample_K_idx(m2),m2)<params.Z_prior_min) || (propos_Z_idx(sample_K_idx(m2),m2)>params.Z_prior_max)
        log_a = -inf;
    else
        propos_B_idx = params.D_prior_B + 0.25*sum([diff(propos_X_idx),diff(propos_Y_idx),diff(propos_Z_idx)].^2 ./ params.dt,[1 2]);
                
        log_a = A*log(sample_B_idx/propos_B_idx);
    
        if ~get_sanity_check(log_a)
            keyboard
        end
    end % log_a
        
    if log_a>=0 || rand < exp(log_a)
        sample_B_idx = propos_B_idx;
        sample_X_idx = propos_X_idx;
        sample_Y_idx = propos_Y_idx;
        sample_Z_idx = propos_Z_idx;
    
        rec(1) = rec(1) + 1;
    end
        rec(2) = rec(2) + 1;

end % rep


%% transfer out
sample_X(:,m_idx) = sample_X_idx;
sample_Y(:,m_idx) = sample_Y_idx;
sample_Z(:,m_idx) = sample_Z_idx;

[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,b,params);

end % function
