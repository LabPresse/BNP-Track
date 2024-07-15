function sample = chainer_init_sample(params)

    %% counter
    sample.i = 0;

    %% temperature for simulated annealing
    sample.T = params.T_init;

    %% diffusion coefficient
    sample.D = params.D_prior_phi * params.D_prior_chi / randg(params.D_prior_phi);

    %% emission rate
    sample.h = params.h_prior_psi / params.h_prior_phi * randg(params.h_prior_phi);

    %% background photon flux
    sample.F = params.F_prior_psi / params.F_prior_phi * randg(params.F_prior_phi);

    %% EMCCD gain
    sample.G = params.G_prior_phi * params.G_prior_chi / randg(params.G_prior_phi);

    %% emitter loads
    sample.bm = (rand(1, params.M) <= params.bm_prior_gamma / (params.bm_prior_gamma + params.M - 1));

    sample.sm = ones(params.N, params.M);

    [sample.Xm, sample.Ym, sample.Zm] = sampler_update_locs_cm(params.M, sample.D, params);

    %% record acceptance ratios for some adapative tuning
    sample.rec_Fh = repmat([0; eps], 1, 3); % emission rates
    sample.rec_XYZ = repmat([0; eps], 1, 2); % locator proposals
end
