function [sample_C,sample_h,sample_K,sample_X,sample_Y,sample_Z,rec] = sampler_update_PPC( ...
          sample_C,sample_h,sample_K,sample_X,sample_Y,sample_Z,rec, ...
          b,F,params)

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

Ap = diff(params.x_bnd)*diff(params.y_bnd);

%%
for n = randperm(params.N)

    if n > 1
        MC_prev =     sample_C(n-1);
        Mh_prev =     sample_h(n-1);
        Mx_prev = sample_X_idx(n-1,:);
        My_prev = sample_Y_idx(n-1,:);
        Mz_prev = sample_Z_idx(n-1,:);
        dt_prev = params.dt(n-1);
    end
    if n < params.N
        MC_next =     sample_C(n+1);
        Mh_next =     sample_h(n+1);
        Mx_next = sample_X_idx(n+1,:);
        My_next = sample_Y_idx(n+1,:);
        Mz_next = sample_Z_idx(n+1,:);
        dt_next = params.dt(n);
    end

    sample_L = get_log_L(sample_C(n),sample_h(n),sample_X_idx(n,:),sample_Y_idx(n,:),sample_Z_idx(n,:),M_idx,[],n,Ap,F,params);

    for rep = 1:poissrnd( 1*M_idx )
    
        % update anchors
        for m = 1:M_idx
            k = find( ( sample_X_idx(:,m)>params.X_prior_min) & ( sample_X_idx(:,m)<params.X_prior_max) ...
                    & ( sample_Y_idx(:,m)>params.Y_prior_min) & ( sample_Y_idx(:,m)<params.Y_prior_max) ...
                    & ( sample_Z_idx(:,m)>params.Z_prior_min) & ( sample_Z_idx(:,m)<params.Z_prior_max) ) ;
            sample_K_idx(m) = k(randi(length(k)));
        end

        idx = find( sample_K_idx==n );
    
        Dr = params.D_prior_B + 0.25*sum( [diff(sample_X_idx),diff(sample_Y_idx),diff(sample_Z_idx)].^2 ./ params.dt,[1 2] );
        Dr = Dr/randg(params.D_prior_A+1.5*(params.N-1)*M_idx);

        tempor_dt    = params.dt;
        tempor_C     = sample_C    ;
        tempor_h     = sample_h    ;
        tempor_X_idx = sample_X_idx;
        tempor_Y_idx = sample_Y_idx;
        tempor_Z_idx = sample_Z_idx;

        if n==1
            tempor_dt(1)    = [];
        elseif n==params.N
            tempor_dt(end)  = [];
        else
            tempor_dt(n-1)  = tempor_dt(n-1) + tempor_dt(n);
            tempor_dt(n) = [];
        end
        tempor_C(n)       = [];
        tempor_h(n)       = [];
        tempor_X_idx(n,:) = [];
        tempor_Y_idx(n,:) = [];
        tempor_Z_idx(n,:) = [];

        dC =                    0.25*sum(  diff(tempor_C    )                                       .^2 ./ tempor_dt,[1 2] );
        dh =                    0.25*sum(  diff(tempor_h    )                                       .^2 ./ tempor_dt,[1 2] );
        dr = params.D_prior_B + 0.25*sum( [diff(tempor_X_idx),diff(tempor_Y_idx),diff(tempor_Z_idx)].^2 ./ tempor_dt,[1 2] );
        
        dC = dC/randg(                 0.5*(params.N-2)      );
        dh = dh/randg(                 0.5*(params.N-2)      );
        dr = dr/randg(params.D_prior_A+1.5*(params.N-2)*M_idx);

        if n==1
            MC = MC_next;
            Mh = Mh_next;
            Mx = Mx_next;
            My = My_next;
            Mz = Mz_next;
            vC = 2*dC*dt_next;
            vh = 2*dh*dt_next;
            vr = 2*dr*dt_next;
            Vr = 2*Dr*dt_next;
        elseif n==params.N
            MC = MC_prev;
            Mh = Mh_prev;
            Mx = Mx_prev;
            My = My_prev;
            Mz = Mz_prev;
            vC = 2*dC*dt_prev;
            vh = 2*dh*dt_prev;
            vr = 2*dr*dt_prev;
            Vr = 2*Dr*dt_prev;
        else
            MC = ( MC_prev*dt_next + MC_next*dt_prev )/(  dt_next + dt_prev );
            Mh = ( Mh_prev*dt_next + Mh_next*dt_prev )/(  dt_next + dt_prev );
            Mx = ( Mx_prev*dt_next + Mx_next*dt_prev )/(  dt_next + dt_prev );
            My = ( My_prev*dt_next + My_next*dt_prev )/(  dt_next + dt_prev );
            Mz = ( Mz_prev*dt_next + Mz_next*dt_prev )/(  dt_next + dt_prev );
            vC =  2*dC/( 1/dt_prev+1/dt_next );
            vh =  2*dh/( 1/dt_prev+1/dt_next );
            vr =  2*dr/( 1/dt_prev+1/dt_next );
            Vr =  2*Dr/( 1/dt_prev+1/dt_next );
        end            

        % pick slice
        log_U = log(rand);
        
        % pick ellipse
        UC = sqrt(vC)*randn(1,1    );
        Uh = sqrt(vh)*randn(1,1    );
        Ux = sqrt(vr)*randn(1,M_idx);
        Uy = sqrt(vr)*randn(1,M_idx);
        Uz = sqrt(vr)*randn(1,M_idx);

        % pick interval
        T_min =       - 2*pi*rand;
        T_max = T_min + 2*pi;


        % keep resampling
        while true 
        
            rec(2) = rec(2) + 1;
            
            % get proposal
            propos_T = T_min + (T_max-T_min)*rand;

            propos_c = MC+(sample_C(n)      -MC)*cos(propos_T)+UC*sin(propos_T);
            propos_H = Mh+(sample_h(n)      -Mh)*cos(propos_T)+Uh*sin(propos_T);
            propos_x = Mx+(sample_X_idx(n,:)-Mx)*cos(propos_T)+Ux*sin(propos_T);
            propos_y = My+(sample_Y_idx(n,:)-My)*cos(propos_T)+Uy*sin(propos_T);
            propos_z = Mz+(sample_Z_idx(n,:)-Mz)*cos(propos_T)+Uz*sin(propos_T);

            propos_L = get_log_L(abs(propos_c),abs(propos_H),propos_x,propos_y,propos_z,M_idx,idx,n,Ap,F,params);

            log_a = propos_L - sample_L ...
                  + 0.5*( (propos_c-MC)^2 - (sample_C(n)-MC)^2 )/vC ...
                  + 0.5*( (propos_H-Mh)^2 - (sample_h(n)-Mh)^2 )/vh ...
                  + 0.5*sum( (propos_x-Mx).^2 - (sample_X_idx(n,:)-Mx).^2 ...
                           + (propos_y-My).^2 - (sample_Y_idx(n,:)-My).^2 ...
                           + (propos_z-Mz).^2 - (sample_Z_idx(n,:)-Mz).^2 )*(1/vr-1/Vr);
            if ~get_sanity_check(log_a)
                keyboard
            end

            % take acceptance test
            if log_U < log_a
                sample_L          =     propos_L;
                sample_C(n)       = abs(propos_c);
                sample_h(n)       = abs(propos_H);
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
function log_L = get_log_L(Cn,hn,X,Y,Z,M,idx,n,Ap,F,params)

if ~isempty(idx) && ( any(X(idx)<params.X_prior_min) || any(X(idx)>params.X_prior_max) ...
                   || any(Y(idx)<params.Y_prior_min) || any(Y(idx)>params.Y_prior_max) ...
                   || any(Z(idx)<params.Z_prior_min) || any(Z(idx)>params.Z_prior_max) ) 
    log_L = -inf;
else
    fn = hn*params.t_exp(n);
    u = Cn*params.t_exp(n)*Ap;
    for m = 1:M
        u = u + get_hG_cnt(fn,X(m),Y(m),Z(m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref); 
    end % m
    log_L = get_log_like(u,F,params,n) ...
          + (params.C_prior_A(n)-1)*log(Cn/params.C_prior_REF(n)) - params.C_prior_A(n)*Cn/params.C_prior_REF(n) ...
          + (params.h_prior_A(n)-1)*log(hn/params.h_prior_REF(n)) - params.h_prior_A(n)*hn/params.h_prior_REF(n);
end

end
