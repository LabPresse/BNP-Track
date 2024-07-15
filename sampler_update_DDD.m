function [K,X,Y,Z,D] = sampler_update_DDD(K,X,Y,Z,b,params)


%% D|active
m = find(b);

D = params.D_prior_B + 0.25*sum( [diff(X(:,m)),diff(Y(:,m)),diff(Z(:,m))].^2 ./ params.dt,[1 2] );
D = D/randg(params.D_prior_A+1.5*(params.N-1)*length(m));


%% inactive|D
m = find(~b);
M = length(m);

if M>0

    K(m) = randi(params.N,1,M);
    k_idx = sub2ind([params.N params.M],K(m),m);
    
    sqrt_two_D_dt = sqrt(2*D*params.dt);

    X(:,m) = cumsum( [zeros(1,M) ; sqrt_two_D_dt.*randn(params.N-1,M)] );
    Y(:,m) = cumsum( [zeros(1,M) ; sqrt_two_D_dt.*randn(params.N-1,M)] );
    Z(:,m) = cumsum( [zeros(1,M) ; sqrt_two_D_dt.*randn(params.N-1,M)] );
    
    X(:,m) = params.X_prior_min + (params.X_prior_max-params.X_prior_min)*rand(1,M) - X(k_idx) + X(:,m);
    Y(:,m) = params.Y_prior_min + (params.Y_prior_max-params.Y_prior_min)*rand(1,M) - Y(k_idx) + Y(:,m);
    Z(:,m) = params.Z_prior_min + (params.Z_prior_max-params.Z_prior_min)*rand(1,M) - Z(k_idx) + Z(:,m);
end

