function [sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_ZZZ( ...
          sample_K,sample_X,sample_Y,sample_Z,rec,...
          b,params)

m_idx = find(b);
M_idx = length(m_idx);

if M_idx<1
    return
end

%% transfer in
sample_X_idx = sample_X(:,m_idx);
sample_Y_idx = sample_Y(:,m_idx);
sample_Z_idx = sample_Z(:,m_idx);
sample_K_idx = sample_K(  m_idx);

%% flip
A = params.D_prior_A+1.5*(params.N-1)*M_idx;
sample_B_idx = params.D_prior_B + 0.25*sum([diff(sample_X_idx),diff(sample_Y_idx),diff(sample_Z_idx)].^2 ./ params.dt,[1 2]);

for rep = 1 : poissrnd(2*params.N*M_idx)

    n = randi(params.N);
    m = randi(M_idx);
    
    propos_z_idx  = sample_Z_idx(:,m);
    propos_z_idx(n) = - propos_z_idx(n);

    if (n==sample_K_idx(m)) && ( (propos_z_idx(n)<params.Z_prior_min) || (propos_z_idx(n)>params.Z_prior_max) )
        log_a = -inf;
    else
        propos_B_idx = sample_B_idx + 0.25*sum(( diff(propos_z_idx     ).^2 ...
                                               - diff(sample_Z_idx(:,m)).^2 ) ./ params.dt );
                
        log_a = A*log(sample_B_idx/propos_B_idx);
    
        if ~get_sanity_check(log_a)
            keyboard
        end
    end % log_a
        
    if log_a>=0 || rand < exp(log_a)
        sample_B_idx      = propos_B_idx;
        sample_Z_idx(:,m) = propos_z_idx;
    
        rec(1) = rec(1) + 1;
    end
        rec(2) = rec(2) + 1;

end % rep


%% transfer out
sample_Z(:,m_idx) = sample_Z_idx;

[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,b,params);

end % function
