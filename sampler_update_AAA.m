function [sample_b,sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_AAA( ...
          sample_b,sample_K,sample_X,sample_Y,sample_Z,rec, ...
          C,h,F,params)
      
[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,sample_b,params);


% get images
m_idx = find(sample_b);
f = h.*params.t_exp;
sample_U = reshape(C.*params.t_exp,1,1,params.N) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
for n = 1:params.N
    for m = m_idx
        sample_U(:,:,n) = sample_U(:,:,n) + get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % m
end % n
sample_L = get_log_like(sample_U,F,params,[]);

sample_A = params.D_prior_A+1.5*(params.N-1)*length(m_idx);
sample_B = params.D_prior_B + 0.25*sum([diff(sample_X(:,m_idx)),diff(sample_Y(:,m_idx)),diff(sample_Z(:,m_idx))].^2 ./ params.dt,[1 2]);


for rep = 1 : poissrnd(params.M)
 
    m = randi(params.M);

    sample_u = sample_U;
    sample_l = sample_L;
    if sample_b(m)
        tempor_A = sample_A - 1.5*(params.N-1);
        tempor_B = sample_B - 0.25*sum([diff(sample_X(:,m)),diff(sample_Y(:,m)),diff(sample_Z(:,m))].^2 ./ params.dt,[1 2]);
        for n = 1:params.N
            sample_U(:,:,n) = sample_U(:,:,n) - get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
        end % n
        sample_L = get_log_like(sample_U,F,params,[]);
    else
        tempor_A = sample_A;
        tempor_B = sample_B;
        for n = 1:params.N
            sample_u(:,:,n) = sample_u(:,:,n) + get_hG_cnt(f(n),sample_X(n,m),sample_Y(n,m),sample_Z(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
        end % n
        sample_l = get_log_like(sample_u,F,params,[]);
    end
    sample_P = sample_L; 
    sample_p = sample_l + params.b_prior_log_p1_m_log_p0;
    sample_Q = max(sample_P,sample_p) + log1p(exp(-abs(sample_P-sample_p)));

    tempor_D = tempor_B/randg(tempor_A);
    sqrt_two_DS_dt = sqrt(2*tempor_D*params.dt);
    propos_K = randi(params.N); %sample_K(m); % randi(params.N);
    propos_X = cumsum( [0 ; sqrt_two_DS_dt.*randn(params.N-1,1)] );
    propos_Y = cumsum( [0 ; sqrt_two_DS_dt.*randn(params.N-1,1)] );
    propos_Z = cumsum( [0 ; sqrt_two_DS_dt.*randn(params.N-1,1)] );
    propos_X = params.X_prior_min + (params.X_prior_max-params.X_prior_min)*rand(1) - propos_X(propos_K) + propos_X;
    propos_Y = params.Y_prior_min + (params.Y_prior_max-params.Y_prior_min)*rand(1) - propos_Y(propos_K) + propos_Y;
    propos_Z = params.Z_prior_min + (params.Z_prior_max-params.Z_prior_min)*rand(1) - propos_Z(propos_K) + propos_Z;

    propos_U = sample_U;
    propos_u = sample_U;
    for n = 1:params.N
        propos_u(:,:,n) = propos_u(:,:,n) + get_hG_cnt(f(n),propos_X(n),propos_Y(n),propos_Z(n),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % n
    propos_l = get_log_like(propos_u,F,params,[]);
    propos_L = sample_L;

    propos_P = propos_L; 
    propos_p = propos_l + params.b_prior_log_p1_m_log_p0;
    propos_Q = max(propos_P,propos_p) + log1p(exp(-abs(propos_P-propos_p)));

    log_a = propos_Q-sample_Q;
    if ~get_sanity_check(log_a)
        keyboard
    end

    if log_a>=0 || rand < exp(log_a)
        sample_K(  m) = propos_K;
        sample_X(:,m) = propos_X;
        sample_Y(:,m) = propos_Y;
        sample_Z(:,m) = propos_Z;
        
        sample_L = propos_L;
        sample_l = propos_l;

        sample_U = propos_U;
        sample_u = propos_u;

        sample_P = propos_P;
        sample_p = propos_p;

    end


    log_a = 1/(1+exp(sample_P-sample_p));
    if ~get_sanity_check(log_a)
        keyboard
    end

    flip = sample_b(m) ~= ( rand<log_a );
    if flip
        sample_b(m) = ~sample_b(m);
        rec(1) = rec(1) + 1;
    end
        rec(2) = rec(2) + 1;

    if sample_b(m)
        sample_U = sample_u;
        sample_L = sample_l;
        
        sample_A = tempor_A + 1.5*(params.N-1);
        sample_B = tempor_B + 0.25*sum([diff(sample_X(:,m)),diff(sample_Y(:,m)),diff(sample_Z(:,m))].^2 ./ params.dt,[1 2]);
    else
        sample_A = tempor_A;
        sample_B = tempor_B;
    end

        
end % rep


[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,sample_b,params);

