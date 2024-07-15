function L = get_log_probs(F,C,b,h,K,X,Y,Z,params)
% L = log( post,like,priors )

if nargin<1
    L = 8;
    return
end

if any(h<0) || any(C<0) 
    error('Non-negative priors are violated')
end

k_idx = sub2ind([params.N params.M],K,1:params.M);
if any( X(k_idx)<params.X_prior_min ) || any( X(k_idx)>params.X_prior_max ) || ...
   any( Y(k_idx)<params.Y_prior_min ) || any( Y(k_idx)>params.Y_prior_max ) || ...
   any( Z(k_idx)<params.Z_prior_min ) || any( Z(k_idx)>params.Z_prior_max )
    error('Position priors are violated')
end

L = zeros(1,get_log_probs);

idx = find(b);
M = length(idx);

%% Likelihood
X_idx = X(:,idx);
Y_idx = Y(:,idx);
Z_idx = Z(:,idx);
f = h.*params.t_exp;
u_cnt = reshape(C.*params.t_exp,1,1,params.N) .* ( diff(params.x_bnd)*diff(params.y_bnd) );
for n = 1:params.N
    for m = 1:M
        u_cnt(:,:,n) = u_cnt(:,:,n) + get_hG_cnt(f(n),X_idx(n,m),Y_idx(n,m),Z_idx(n,m),params.x_bnd,params.y_bnd,params.s_ref,params.z_ref);
    end % m
end % n
L(2) = get_log_like(u_cnt,F,params,[]);


%% Priors
L(3) = M*params.b_prior_log_p1_m_log_p0;
L(4) = sum((params.C_prior_A-1).*log(C./params.C_prior_REF)) - sum(params.C_prior_A.*C./params.C_prior_REF);
L(5) = sum((params.h_prior_A-1).*log(h./params.h_prior_REF)) - sum(params.h_prior_A.*h./params.h_prior_REF);

B = params.D_prior_B + 0.25*sum([diff(X),diff(Y),diff(Z)].^2 ./ params.dt,[1 2]);
L(6) = ( params.D_prior_A + 1.5*(params.N-1)*params.M ) * log(params.D_prior_B/B);


B = params.D_prior_B + 0.25*sum([diff(X_idx),diff(Y_idx),diff(Z_idx)].^2 ./ params.dt,[1 2]);
L(7) = ( params.D_prior_A + 1.5*(params.N-1)*M ) * log(params.D_prior_B/B);

%% Posterior
L(1) = L(2) + L(3) + L(4) + L(5) + L(6); % full with D marginalized
L(8) = L(2)        + L(4) + L(5) + L(7); % conditional on b with D and inactive marginalized


