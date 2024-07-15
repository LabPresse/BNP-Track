function [sample_C,rec] = sampler_update_BCK( ...
          sample_C,h,b,X,Y,Z,F,rec, ...
          params)

m_idx = find(b);
M_idx = length(m_idx);
X_idx = X(:,m_idx);
Y_idx = Y(:,m_idx);
Z_idx = Z(:,m_idx);

for n = 1:params.N

    f = h(n)*params.t_exp(n);
    A_tau_exp = params.t_exp(n) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
    u_cnt = zeros(params.Px,params.Py);
    for m = 1:M_idx
        u_cnt = u_cnt + get_hG_cnt(f,X_idx(n,m),Y_idx(n,m),Z_idx(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);
    end % m

    sample_L = get_L(n,u_cnt,sample_C(n),A_tau_exp,F,params);

    for rep = 1:poissrnd( params.BKC_rep )

        % pick slice
        U_prop = log(rand);

        % get intervals
        c_min = sample_C(n) - params.C_WIN(n) * rand(); 
        c_max =       c_min + params.C_WIN(n);

        i = params.I_max;

                dL_min = get_L(n,u_cnt,abs(c_min),A_tau_exp,F,params) - sample_L;
                dL_max = get_L(n,u_cnt,abs(c_max),A_tau_exp,F,params) - sample_L;
        while i>0 && ( U_prop<dL_min || U_prop<dL_max )
            if rand < 0.5
                c_min = 2*c_min - c_max;
                dL_min = get_L(n,u_cnt,abs(c_min),A_tau_exp,F,params) - sample_L;
            else
                c_max = 2*c_max - c_min;
                dL_max = get_L(n,u_cnt,abs(c_max),A_tau_exp,F,params) - sample_L;
            end
            i = i-1;
            rec(3) = rec(3) + 1;
        end


        while true 

            rec(2) = rec(2) + 1;

            propos_c = c_min + (c_max-c_min)*rand();

            propos_L = get_L(n,u_cnt,abs(propos_c),A_tau_exp,F,params);
            
            log_a = propos_L-sample_L;
            if ~get_sanity_check(log_a)
                keyboard
            end


            % carry acc test
            if U_prop < log_a
                sample_L = propos_L;
                sample_C(n) = abs(propos_c);

                rec(1) = rec(1) + 1;
                break % while true
            else
                % update intervals
                if propos_c < sample_C(n)
                    c_min = propos_c;
                else
                    c_max = propos_c;
                end

            end % acc

        end % while true

    end % rep
end % n


end % function


%% helpers
function log_L = get_L(n,u_cnt,C,A_tau_exp,F,params)

log_L = get_log_like( u_cnt + C*A_tau_exp , F, params , n ) ...
      + (params.C_prior_A(n)-1)*log(C/params.C_prior_REF(n)) - params.C_prior_A(n)*C/params.C_prior_REF(n);

end