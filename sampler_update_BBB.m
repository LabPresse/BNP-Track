function [sample_b,sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_BBB( ...
          sample_b,sample_K,sample_X,sample_Y,sample_Z,rec, ...
          C,h,F,params)
      

% get images
m_idx = find(sample_b);
f = h.*params.t_exp;
u_cnt = reshape(C.*params.t_exp,1,1,params.N) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
for n = 1:params.N
    for m = m_idx
        u_cnt(:,:,n) = u_cnt(:,:,n) + get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % m
end % n
l_cnt = get_log_like(u_cnt,F,params,[]);



for rep = 1 : poissrnd(params.M)

    m = randi(params.M);

    if ~sample_b(m)
        m_idx = find(sample_b);

        D = params.D_prior_B + 0.25*sum( [diff(sample_X(:,m_idx)),diff(sample_Y(:,m_idx)),diff(sample_Z(:,m_idx))].^2 ./ params.dt,[1 2] );
        D = D/randg(params.D_prior_A+1.5*(params.N-1)*length(m_idx));

        sample_K(m) = randi(params.N);

        sqrt_two_D_dt = sqrt(2*D*params.dt);
    
        sample_X(:,m) = cumsum( [0 ; sqrt_two_D_dt.*randn(params.N-1,1)] );
        sample_Y(:,m) = cumsum( [0 ; sqrt_two_D_dt.*randn(params.N-1,1)] );
        sample_Z(:,m) = cumsum( [0 ; sqrt_two_D_dt.*randn(params.N-1,1)] );
        
        sample_X(:,m) = params.X_prior_min + (params.X_prior_max-params.X_prior_min)*rand - sample_X(sample_K(m),m) + sample_X(:,m);
        sample_Y(:,m) = params.Y_prior_min + (params.Y_prior_max-params.Y_prior_min)*rand - sample_Y(sample_K(m),m) + sample_Y(:,m);
        sample_Z(:,m) = params.Z_prior_min + (params.Z_prior_max-params.Z_prior_min)*rand - sample_Z(sample_K(m),m) + sample_Z(:,m);
    end

    U_cnt = u_cnt;
    L_cnt = l_cnt;

    if sample_b(m)
        % remove
        for n = 1:params.N
            U_cnt(:,:,n) = U_cnt(:,:,n) - get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
        end % n
        L_cnt = get_log_like(U_cnt,F,params,[]);
    else
        % add
        for n = 1:params.N
            u_cnt(:,:,n) = u_cnt(:,:,n) + get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);
        end % n
        l_cnt = get_log_like(u_cnt,F,params,[]);
    end % if
    
    log_a = 1/(1+exp( L_cnt-l_cnt-params.b_prior_log_p1_m_log_p0 ));
    if ~get_sanity_check(log_a)
        keyboard
    end

    if sample_b(m) ~= (rand<log_a)
        sample_b(m) = ~sample_b(m);
        rec(1) = rec(1)+1;
    end
        rec(2) = rec(2)+1;

    if ~sample_b(m)
        u_cnt = U_cnt;
        l_cnt = L_cnt;
    end
        
end % rep

[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,sample_b,params);

