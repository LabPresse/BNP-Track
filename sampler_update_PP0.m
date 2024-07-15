function [sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_PP0( ...
          sample_K,sample_X,sample_Y,sample_Z,rec, ...
          b,C,h,F,params)

m_idx = find(b);

if isempty(m_idx)
    return
end

M_idx = length(m_idx);

%% transfer in
sample_X_idx = sample_X(:,m_idx);
sample_Y_idx = sample_Y(:,m_idx);
sample_Z_idx = sample_Z(:,m_idx);
sample_K_idx = sample_K(  m_idx);


%% images and likelihoods
sample_L = nan(params.N,1);

f = h.*params.t_exp;
sample_u_cnt = reshape(C.*params.t_exp,1,1,params.N) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
for n = 1:params.N
    for m = 1:M_idx
        sample_u_cnt(:,:,n) = sample_u_cnt(:,:,n) + get_hG_cnt(f(n),sample_X_idx(n,m),sample_Y_idx(n,m),sample_Z_idx(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % m
    sample_L(n) = get_log_like(sample_u_cnt(:,:,n),F,params,n);
end % n

A = params.D_prior_A+1.5*(params.N-1)*M_idx;

sample_B = params.D_prior_B + 0.25*sum([diff(sample_X_idx),diff(sample_Y_idx),diff(sample_Z_idx)].^2 ./ params.dt,[1 2]);


%%
for rep = 1:poissrnd(5*M_idx*params.N*params.PP0_rep)

    i = randi(5);
    m = randi(M_idx);
    n = randi(params.N);


    k = find( ( sample_X_idx(:,m)>params.X_prior_min) & ( sample_X_idx(:,m)<params.X_prior_max) ...
            & ( sample_Y_idx(:,m)>params.Y_prior_min) & ( sample_Y_idx(:,m)<params.Y_prior_max) ...
            & ( sample_Z_idx(:,m)>params.Z_prior_min) & ( sample_Z_idx(:,m)<params.Z_prior_max) ) ;
    sample_K_idx(m) = k(randi(length(k)));
    
    [sample_X_idx(:,m),sample_Y_idx(:,m),sample_Z_idx(:,m),sample_L(n),sample_u_cnt(:,:,n),sample_B,rec] = update_single( ...
     sample_X_idx(:,m),sample_Y_idx(:,m),sample_Z_idx(:,m),sample_L(n),sample_u_cnt(:,:,n),sample_B,rec, ...
     sample_K_idx(  m),f(n),i,n,A,F,params);

    k = find( ( sample_X_idx(:,m)>params.X_prior_min) & ( sample_X_idx(:,m)<params.X_prior_max) ...
            & ( sample_Y_idx(:,m)>params.Y_prior_min) & ( sample_Y_idx(:,m)<params.Y_prior_max) ...
            & ( sample_Z_idx(:,m)>params.Z_prior_min) & ( sample_Z_idx(:,m)<params.Z_prior_max) ) ;
    sample_K_idx(m) = k(randi(length(k)));

end % rep


%% transfer out
sample_X(:,m_idx) = sample_X_idx;
sample_Y(:,m_idx) = sample_Y_idx;
sample_Z(:,m_idx) = sample_Z_idx;
sample_K(  m_idx) = sample_K_idx;

[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,b,params);



end % function



%%
function [sample_X,sample_Y,sample_Z,sample_L,sample_u,sample_B,rec] = update_single( ...
          sample_X,sample_Y,sample_Z,sample_L,sample_u,sample_B,rec,...
          K,f,i,n,A,F,params)

% remove
sample_u = sample_u ...
         - get_hG_cnt(f,sample_X(n),sample_Y(n),sample_Z(n),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);

% pick slice
log_U = log(rand);

% choose brackets
if n~=K
    switch i
        case 1
            X_win = params.x_win;
            Y_win = params.y_win;
            Z_win = params.z_win;
        case 2
            X_win = params.x_win;
            Y_win = params.y_win;
            Z_win = 0;
        case 3
            X_win = 0;
            Y_win = 0;
            Z_win = params.z_win;
        case 4
            X_win = params.x_win;
            Y_win = 0;
            Z_win = 0;
        case 5
            X_win = 0;
            Y_win = params.y_win;
            Z_win = 0;
    end % switch
    x_min = sample_X(n) - X_win * rand(); 
    y_min = sample_Y(n) - Y_win * rand(); 
    z_min = sample_Z(n) - Z_win * rand(); 
    x_max = x_min + X_win;
    y_max = y_min + Y_win;
    z_max = z_min + Z_win;
else
    switch i
        case 1
            x_min = params.X_prior_min;
            y_min = params.Y_prior_min;
            z_min = params.Z_prior_min;
            x_max = params.X_prior_max;
            y_max = params.Y_prior_max;
            z_max = params.Z_prior_max;
        case 2
            x_min = params.X_prior_min;
            y_min = params.Y_prior_min;
            z_min = sample_Z(n);
            x_max = params.X_prior_max;
            y_max = params.Y_prior_max;
            z_max = sample_Z(n);
        case 3
            x_min = sample_X(n);
            y_min = sample_Y(n);
            z_min = params.Z_prior_min;
            x_max = sample_X(n);
            y_max = sample_Y(n);
            z_max = params.Z_prior_max;
        case 4
            x_min = params.X_prior_min;
            y_min = sample_Y(n);
            z_min = sample_Z(n);
            x_max = params.X_prior_max;
            y_max = sample_Y(n);
            z_max = sample_Z(n);
        case 5
            x_min = sample_X(n);
            y_min = params.Y_prior_min;
            z_min = sample_Z(n);
            x_max = sample_X(n);
            y_max = params.Y_prior_max;
            z_max = sample_Z(n);
    end % switch
end % if


propos_X = sample_X;
propos_Y = sample_Y;
propos_Z = sample_Z;        

% keep resampling
while true 

    rec(2) = rec(2) + 1;
    
    % get proposal
    propos_X(n) = x_min + (x_max-x_min)*rand();
    propos_Y(n) = y_min + (y_max-y_min)*rand();
    propos_Z(n) = z_min + (z_max-z_min)*rand();

    % add
    propos_u = sample_u ...
             + get_hG_cnt( f , ...
                           propos_X(n),propos_Y(n),propos_Z(n),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);
    propos_L = get_log_like(propos_u,F,params,n);
    propos_B = sample_B + 0.25*sum( ( [diff(propos_X),diff(propos_Y),diff(propos_Z)].^2 ...
                                    - [diff(sample_X),diff(sample_Y),diff(sample_Z)].^2 )./params.dt,[1 2]);
    log_a = propos_L - sample_L ...
          + A*log(sample_B/propos_B);
    if ~get_sanity_check(log_a)
        keyboard
    end

    % take acceptance test
    if log_U < log_a
        sample_u = propos_u;
        sample_L = propos_L;
        sample_X = propos_X;
        sample_Y = propos_Y;
        sample_Z = propos_Z;
        sample_B = propos_B;
        
        rec(1) = rec(1) + 1;
        break % while true
    else
        if propos_X(n) < sample_X(n)
            x_min = propos_X(n);
        else
            x_max = propos_X(n);
        end
        if propos_Y(n) < sample_Y(n)
            y_min = propos_Y(n);
        else
            y_max = propos_Y(n);
        end
        if propos_Z(n) < sample_Z(n)
            z_min = propos_Z(n);
        else
            z_max = propos_Z(n);
        end
    end % acc
    
end % while

end % function
