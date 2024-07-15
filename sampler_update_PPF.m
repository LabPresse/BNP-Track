function [sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_PPF( ...
          sample_K,sample_X,sample_Y,sample_Z,rec, ...
          C,h,b,F,params)

m_idx = find(b);

if isempty(m_idx)
    return
end

M_idx = length(m_idx);


%% transfer in
sample_K_idx = sample_K(  m_idx);
sample_X_idx = sample_X(:,m_idx);
sample_Y_idx = sample_Y(:,m_idx);
sample_Z_idx = sample_Z(:,m_idx);


%% images and likelihoods

for n = randperm(params.N)

    Cn_t_exp_Ap = C(n)*params.t_exp(n) * ( diff(params.x_bnd)*diff(params.y_bnd) );
    fn          = h(n)*params.t_exp(n) ;

    if n > 1
        Mx_prev = sample_X_idx(n-1,:);
        My_prev = sample_Y_idx(n-1,:);
        Mz_prev = sample_Z_idx(n-1,:);
        dt_prev = params.dt(n-1);
    end
    if n < params.N
        Mx_next = sample_X_idx(n+1,:);
        My_next = sample_Y_idx(n+1,:);
        Mz_next = sample_Z_idx(n+1,:);
        dt_next = params.dt(n);
    end

    sample_L = get_log_L(fn,sample_X_idx(n,:),sample_Y_idx(n,:),sample_Z_idx(n,:),M_idx,[],n,Cn_t_exp_Ap,F,params);

    for rep = 1:poissrnd( M_idx*params.PPF_rep )
    
        % update anchors
        for m = 1:M_idx
            k = find( ( sample_X_idx(:,m)>params.X_prior_min) & ( sample_X_idx(:,m)<params.X_prior_max) ...
                    & ( sample_Y_idx(:,m)>params.Y_prior_min) & ( sample_Y_idx(:,m)<params.Y_prior_max) ...
                    & ( sample_Z_idx(:,m)>params.Z_prior_min) & ( sample_Z_idx(:,m)<params.Z_prior_max) ) ;
            sample_K_idx(m) = k(randi(length(k)));
        end

        idx = find( sample_K_idx==n );
    
        Dr = params.D_prior_B + get_S_without([],[sample_X_idx,sample_Y_idx,sample_Z_idx],params.dt);
        Dr = Dr/randg(params.D_prior_A+1.5*(params.N-1)*M_idx);

        dr = params.D_prior_B + get_S_without(n,[sample_X_idx,sample_Y_idx,sample_Z_idx],params.dt);
        dr = dr/randg(params.D_prior_A+1.5*(params.N-2)*M_idx);

        if n==1
            Mx = Mx_next;
            My = My_next;
            Mz = Mz_next;
            Vr = 2*Dr*dt_next;
            vr = 2*dr*dt_next;
        elseif n==params.N
            Mx = Mx_prev;
            My = My_prev;
            Mz = Mz_prev;
            Vr = 2*Dr*dt_prev;
            vr = 2*dr*dt_prev;
        else
            Mx = ( Mx_prev*dt_next + Mx_next*dt_prev )/(  dt_next + dt_prev );
            My = ( My_prev*dt_next + My_next*dt_prev )/(  dt_next + dt_prev );
            Mz = ( Mz_prev*dt_next + Mz_next*dt_prev )/(  dt_next + dt_prev );
            Vr = 2*Dr/( 1/dt_prev+1/dt_next );
            vr = 2*dr/( 1/dt_prev+1/dt_next );
        end            

        % pick slice
        log_U = log(rand);

        % pick ellipse
        ux = sqrt(vr)*randn(1,M_idx);
        uy = sqrt(vr)*randn(1,M_idx);
        uz = sqrt(vr)*randn(1,M_idx);

        % pick interval
        T_min =       - 2*pi*rand;
        T_max = T_min + 2*pi;


        % keep resampling
        while true 
        
            rec(2) = rec(2) + 1;
            
            % get proposal
            propos_T = T_min + (T_max-T_min)*rand;

            propos_x = Mx+(sample_X_idx(n,:)-Mx)*cos(propos_T)+ux*sin(propos_T);
            propos_y = My+(sample_Y_idx(n,:)-My)*cos(propos_T)+uy*sin(propos_T);
            propos_z = Mz+(sample_Z_idx(n,:)-Mz)*cos(propos_T)+uz*sin(propos_T);

            propos_L = get_log_L(fn,propos_x,propos_y,propos_z,M_idx,idx,n,Cn_t_exp_Ap,F,params);

            log_a = propos_L - sample_L ...
                  + 0.5*sum( (propos_x-Mx).^2 - (sample_X_idx(n,:)-Mx).^2 ...
                           + (propos_y-My).^2 - (sample_Y_idx(n,:)-My).^2 ...
                           + (propos_z-Mz).^2 - (sample_Z_idx(n,:)-Mz).^2 )*(1/vr-1/Vr);
            if ~get_sanity_check(log_a)
                keyboard
            end

            % take acceptance test
            if log_U < log_a
                sample_L          =     propos_L;
                sample_X_idx(n,:) =     propos_x;
                sample_Y_idx(n,:) =     propos_y;
                sample_Z_idx(n,:) =     propos_z;
                
                rec(1) = rec(1) + 1;
                break % while true
            else
                if propos_T<0
                    T_min = propos_T;
                else
                    T_max = propos_T; 
                end
            end % acc
            
        end % while

    end % rep

end % n

% re-update anchors
for m = 1:M_idx
    k = find( ( sample_X_idx(:,m)>params.X_prior_min) & ( sample_X_idx(:,m)<params.X_prior_max) ...
            & ( sample_Y_idx(:,m)>params.Y_prior_min) & ( sample_Y_idx(:,m)<params.Y_prior_max) ...
            & ( sample_Z_idx(:,m)>params.Z_prior_min) & ( sample_Z_idx(:,m)<params.Z_prior_max) ) ;
    sample_K_idx(m) = k(randi(length(k)));
end

%% transfer out
sample_K(  m_idx) = sample_K_idx;
sample_X(:,m_idx) = sample_X_idx;
sample_Y(:,m_idx) = sample_Y_idx;
sample_Z(:,m_idx) = sample_Z_idx;

[sample_K,sample_X,sample_Y,sample_Z] = sampler_update_DDD(sample_K,sample_X,sample_Y,sample_Z,b,params);


end % function


%% ------------------------------------------------------------------------
function log_L = get_log_L(fn,X,Y,Z,M,idx,n,Cn_t_exp_Ap,F,params)

if ~isempty(idx) && ( any(X(idx)<params.X_prior_min) || any(X(idx)>params.X_prior_max) ...
                   || any(Y(idx)<params.Y_prior_min) || any(Y(idx)>params.Y_prior_max) ...
                   || any(Z(idx)<params.Z_prior_min) || any(Z(idx)>params.Z_prior_max) ) 
    log_L = -inf;
else
    u = Cn_t_exp_Ap;
    for m = 1:M
        u = u + get_hG_cnt(fn,X(m),Y(m),Z(m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % m
    log_L = get_log_like(u,F,params,n);
end

end


%%
function S = get_S_without(n,f,dt)

if ~isempty(n)
    t = cumsum([0;dt]);
    f(n,:) = [];
    t(n  ) = [];
    dt = diff(t);
end

S = 0.25*sum( diff(f).^2 ./ dt , 'all' );

end
