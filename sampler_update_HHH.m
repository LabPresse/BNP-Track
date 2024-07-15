function [sample_C,sample_h,rec] = sampler_update_HHH( ...
          sample_C,sample_h,rec,b,X,Y,Z,F, ...
          params)

m_idx = find(b);
if isempty(m_idx)
    sample_h = params.h_prior_REF./params.h_prior_A.*randg(params.h_prior_A,params.N,1);
    return
end

M_idx = length(m_idx);
X_idx = X(:,m_idx);
Y_idx = Y(:,m_idx);
Z_idx = Z(:,m_idx);

for n = 1:params.N

    V_cnt = params.t_exp(n) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
    U_cnt = zeros(params.Px,params.Py);
    for m = 1:M_idx
        U_cnt = U_cnt + get_hG_cnt(params.t_exp(n),X_idx(n,m),Y_idx(n,m),Z_idx(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);
    end % m

    sample_L = get_L(sample_C(n),sample_h(n),n,V_cnt,U_cnt,F,params);

    for rep = 1:poissrnd( params.HHH_rep )

        % pick slice
        U_prop = log(rand);

        % get intervals
        c_min = sample_C(n) - params.C_WIN(n) * rand; 
        c_max =       c_min + params.C_WIN(n);
        
        H_min = sample_h(n) - params.h_WIN(n) * rand; 
        H_max =       H_min + params.h_WIN(n);

        % adapt intervals
        i = params.I_max;

                dL_min_min = get_L(abs(c_min),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_min_max = get_L(abs(c_min),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_max_min = get_L(abs(c_max),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_max_max = get_L(abs(c_max),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
        while i>0 && U_prop < max([dL_min_min,dL_min_max,dL_max_min,dL_max_max])
            t = rand;
            if     t < 0.25
                c_min = 2*c_min - c_max;
                dL_min_min = get_L(abs(c_min),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_min_max = get_L(abs(c_min),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
            elseif t < 0.50
                c_max = 2*c_max - c_min;
                dL_max_min = get_L(abs(c_max),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_max_max = get_L(abs(c_max),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
            elseif t < 0.75
                H_min = 2*H_min - H_max;
                dL_min_min = get_L(abs(c_min),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_max_min = get_L(abs(c_max),abs(H_min),n,V_cnt,U_cnt,F,params) - sample_L;
            else
                H_max = 2*H_max - H_min;
                dL_min_max = get_L(abs(c_min),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
                dL_max_max = get_L(abs(c_max),abs(H_max),n,V_cnt,U_cnt,F,params) - sample_L;
            end
            i = i-1;
            rec(3) = rec(3) + 1;
        end


        while true 

            rec(2) = rec(2) + 1;

            propos_c = c_min + (c_max-c_min)*rand;
            propos_H = H_min + (H_max-H_min)*rand;

            propos_L = get_L(abs(propos_c),abs(propos_H),n,V_cnt,U_cnt,F,params);
            
            log_a = propos_L-sample_L;
            if ~get_sanity_check(log_a)
                keyboard
            end


            % carry acc test
            if U_prop < log_a
                sample_L    =     propos_L ;
                sample_C(n) = abs(propos_c);
                sample_h(n) = abs(propos_H);

                rec(1) = rec(1) + 1;
                break % while true
            else
                % update intervals
                if propos_c < sample_C(n)
                    c_min = propos_c;
                else
                    c_max = propos_c;
                end
                if propos_H < sample_h(n)
                    H_min = propos_H;
                else
                    H_max = propos_H;
                end
            end % acc

        end % while true

    end % rep
end % n


end % function


%% helpers

function log_L = get_L(C,h,n,V_cnt,U_cnt,F,params)

log_L = get_log_like( C*V_cnt+h*U_cnt ,F,params,n) ...
      + (params.h_prior_A(n)-1)*log(h/params.h_prior_REF(n)) - params.h_prior_A(n)*h/params.h_prior_REF(n) ...
      + (params.C_prior_A(n)-1)*log(C/params.C_prior_REF(n)) - params.C_prior_A(n)*C/params.C_prior_REF(n);

end

