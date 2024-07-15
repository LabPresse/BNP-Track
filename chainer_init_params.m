function params = chainer_init_params(opts)


%% Samplers
params.F_IC = size(opts.w_cnt,1)*size(opts.w_cnt,2);
params.F_TC = 1;
params.i_sc = 10/log(2);


%% Set-up
params.units = opts.units;

params.t_ded = double( opts.t_ded );  % [time] dead time 

params.x_bnd = double( opts.x_bnd );  % [length] pixel edges 
params.y_bnd = double( opts.y_bnd );  % [length] pixel edges
params.t_exp = double( opts.t_exp );  % [time] exposure times

params.wM = double( opts.wM ); % [image]          read-out offset
params.wV = double( opts.wV ); % [image]^2        read-out variance
params.wG = double( opts.wG ); % [image]/[photon] overall gain
params.wF = double( opts.wF ); %                  excess noise factor

params.Px = size(opts.w_cnt,1);
params.Py = size(opts.w_cnt,2);
params.N  = size(opts.w_cnt,3);
params.M  = 25;

params.s_ref = double( 0.21*opts.lambda/opts.nNA );                 % [length]
params.z_ref = double( 4*pi*opts.nRI*params.s_ref^2/opts.lambda );  % [length]

%% Consistency check
params.x_bnd = reshape(params.x_bnd,params.Px+1,1);
params.y_bnd = reshape(params.y_bnd,1,params.Py+1);
params.t_exp = reshape(params.t_exp,params.N,1   );


%% Pre-processed measurments
params.dW_cnt = ( double( opts.w_cnt ) - params.wM )/params.wG; % [photons]

if any(diff(params.x_bnd)<0) || ...
   any(diff(params.y_bnd)<0)
    error('Inconsistent arrangement. Sort x_bnd, y_bnd')
end

params.t_bnd = [0;cumsum(params.t_exp+params.t_ded)];
params.t_mid = 0.5*(params.t_bnd(1:end-1)+params.t_bnd(2:end));

params.dt = diff(params.t_mid);


%% Dynamics
D_ref = 0.01; % [area]/[time];
params.D_prior_A = 2;
params.D_prior_B = (params.D_prior_A-1)*D_ref;


%% Positions
ds = 1; % [length]
dz = ds + 0.5*min(params.x_bnd(end)-params.x_bnd(1),params.y_bnd(end)-params.y_bnd(1));

params.X_prior_min = params.x_bnd(  1) - ds;    % [lenght]
params.Y_prior_min = params.y_bnd(  1) - ds;    % [lenght]
params.Z_prior_min =                   - dz;    % [lenght]

params.X_prior_max = params.x_bnd(end) + ds;    % [lenght]
params.Y_prior_max = params.y_bnd(end) + ds;    % [lenght]
params.Z_prior_max =                   + dz;    % [lenght]

%% Photon emission rates
tot_pht = reshape(sum(params.dW_cnt,[1,2]),params.N,1,1);
tot_are = (params.x_bnd(end)-params.x_bnd(1))*(params.y_bnd(end)-params.y_bnd(1));

params.C_prior_A   = tot_pht;
params.C_prior_REF = tot_pht/tot_are./params.t_exp;

area_ref = 1; % [micron]^2
params.h_prior_A   = repmat(2,params.N,1);
params.h_prior_REF = tot_pht./params.t_exp/tot_are*area_ref;

%% Loads
log_b_prior_gamma = 0;

if params.M <= exp( log_b_prior_gamma )
    error('Inconsistent BNP approximation. Adjust M or \gamma')
end

b_prior_log_p1 = log_b_prior_gamma - log(params.M);
b_prior_log_p0 = log1p(-exp(b_prior_log_p1));
params.b_prior_log_p1_m_log_p0 = b_prior_log_p1 - b_prior_log_p0;

%% ground
if isfield(opts,'ground')
    params.ground = opts.ground;
end

%% MCMC Samplers

params.x_win = params.X_prior_max-params.X_prior_min;
params.y_win = params.Y_prior_max-params.Y_prior_min;
params.z_win = params.Z_prior_max-params.Z_prior_min;

params.I_max = inf;

params.C_WIN = 0.01*params.C_prior_REF;
params.h_WIN = 0.10*params.h_prior_REF;

params.BKC_rep = 25;
params.HHH_rep = 25;
params.PP0_rep =  5;
params.PPF_rep =  1;

