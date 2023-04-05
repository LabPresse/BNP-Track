function params = chainer_init_params(opts)

    %% Units
    params.units = opts.units;

    %% Sizes
    params.Px = size(opts.w_cnt, 1);
    params.Py = size(opts.w_cnt, 2);
    params.P = params.Px * params.Py;
    params.N = size(opts.w_cnt, 3);

    params.M = 100;
    params.K = 5; % points per frame

    %% Pixel boundaries and data
    params.x_bnd = reshape(opts.x_bnd, params.Px + 1, 1); % [length]
    params.y_bnd = reshape(opts.y_bnd, 1, params.Py + 1); % [length]
    params.px_area = diff(params.x_bnd) * diff(params.y_bnd);
    params.px_area_col = reshape(diff(params.x_bnd) * diff(params.y_bnd), params.P, 1);

    params.w_cnt = reshape(opts.w_cnt, params.Px, params.Py, params.N);

    %% Image discretization
    params.ak = [0.5 ones(1, params.K - 2) 0.5] / (params.K - 1); % composite trapezoid

    %% Time levels
    params.dt_stp = opts.dt_stp; % [time] frame separation
    params.dt_exp = opts.dt_exp; % [time] exposure time

    params.t_bnd = params.dt_stp * (0:params.N)';
    params.t_mid = reshape(bsxfun( ...
        @plus, ...
        params.t_bnd(1:params.N), ...
        linspace( ...
        0.5 * (params.dt_stp - params.dt_exp), ...
        params.dt_stp - 0.5 * (params.dt_stp - params.dt_exp), ...
        params.K ...
    ) ...
    )', params.N * params.K, 1);
    params.t_idx = reshape(1:params.N * params.K, params.K, params.N)';

    %% Diffusion coefficient
    params.D_prior_phi = 5;
    params.D_prior_chi = (1 - 1 / params.D_prior_phi) * 0.1; % [area/time]

    %% Optical parameters
    params.PSF_params = opts.PSF_params;

    params.f = opts.f; % excess noise factor
    params.w_cnt_log = reshape(log(params.w_cnt / params.f), params.P, params.N);

    params.F_prior_phi = 2;
    params.F_prior_psi = params.f ...
        * mean(params.w_cnt(:))^2 ...
        / var(params.w_cnt(:)) ...
        / params.dt_exp ...
        / (sum(params.px_area, 'all') / params.Px / params.Py); % photons/[area*time]

    params.h_prior_phi = 2;
    params.h_prior_psi = params.F_prior_psi * 5 * pi * 2 * log(2) * params.PSF_params.s_ref^2; % photons/[time]

    params.G_prior_phi = 2;
    params.G_prior_chi = (1 - 1 / params.G_prior_phi) ...
        * (var(params.w_cnt(:)) / mean(params.w_cnt(:)) / params.f); % [image]/photon

    %% Beta-Bernoulli
    params.bm_prior_gamma = 0.05;

    %% Location priors
    params.Xm_prior_mu = 0.5 * (params.x_bnd(1) + params.x_bnd(end)); % [length]
    params.Ym_prior_mu = 0.5 * (params.y_bnd(1) + params.y_bnd(end)); % [length]
    % params.Zm_prior_mu = params.PSF_params.z_ref; % [length]
    params.Zm_prior_mu = 0;

    params.Xm_prior_sg = 0.3 * (params.x_bnd(end) - params.x_bnd(1)); % [length]
    params.Ym_prior_sg = 0.3 * (params.y_bnd(end) - params.y_bnd(1)); % [length]
    params.Zm_prior_sg = params.PSF_params.z_ref; % [length]

    %% aux sampler's params

    params.i_skip = 1;

    params.T_init = 100;
    %params.T_pace = 0.01;
    params.T_end = 1000;

    params.MH_sc = [[100 100 2] [0.05 0.05 0.05]];
    %[[100 100 2] [0.1 * 3 * params.K 1 0.01]];

    %% misc
    % ground truth
    if isfield(opts, 'ground')
        params.ground = opts.ground;
    end
