function [bkndflux, emisrate, gain, rec, params] = ...
        sampler_update_bknd_and_gain( ...
        bkndflux, ...
        emisrate, ...
        active, ...
        activeX, ...
        activeY, ...
        activeZ, ...
        rec, ...
        T, ...
        params ...
    )

    px_PSF = zeros(params.Px, params.Py, params.N);

    for n = 1:params.N

        for m = 1:nnz(active)

            for k = 1:params.K
                px_PSF(:, :, n) = px_PSF(:, :, n) ...
                    + params.ak(k) * get_g_PSF( ...
                    activeX(params.t_idx(n, k), m), ...
                    activeY(params.t_idx(n, k), m), ...
                    activeZ(params.t_idx(n, k), m), ...
                    params.x_bnd, ...
                    params.y_bnd, ...
                    params.PSF_params);
            end

        end

    end

    effective_cnt = params.w_cnt / params.f;
    total_cnt = sum(effective_cnt, 'all');
    log_eff_cnt = log(effective_cnt);

    %%
    effective_photon = params.dt_exp / params.f * (bkndflux * params.px_area + emisrate * px_PSF);
    total_photon = sum(effective_photon, 'all');
    sample_E = sum(effective_photon .* log_eff_cnt - gammaln(effective_photon), 'all');

    for rep = 1:75

        % draw proposals

        [bkndflux_prop, log_a] = get_proposal(bkndflux, 0, params.MH_sc(1), 1);
        emisrate_prop = emisrate;

        effective_photon_prop = params.dt_exp / params.f * (bkndflux_prop * params.px_area + emisrate_prop * px_PSF);
        total_photon_prop = sum(effective_photon_prop, 'all');
        propos_E = sum(effective_photon_prop .* log_eff_cnt - gammaln(effective_photon_prop), 'all');

        log_a = log_a ...
            + (params.F_prior_phi - 1) * log(bkndflux_prop / bkndflux) ...
            + params.F_prior_phi * (bkndflux - bkndflux_prop) / params.F_prior_psi ...
            + (propos_E - sample_E) / T ...
            + (total_photon - total_photon_prop) /T * log(params.G_prior_chi * params.G_prior_phi + total_cnt / T) ...
            + gammaln(params.G_prior_phi + total_photon_prop / T) ...
            - gammaln(params.G_prior_phi + total_photon / T);

        if log(rand) < log_a
            bkndflux = bkndflux_prop;
            total_photon = total_photon_prop;
            sample_E = propos_E;

            rec(1, 1) = rec(1, 1) + 1;
        end

        rec(2, 1) = rec(2, 1) + 1;

    end

    %% update G
    gain = (params.G_prior_chi * params.G_prior_phi + total_cnt / T) ...
        / randg(params.G_prior_phi + total_photon / T);

end

function [x_prop, log_a] = get_proposal(x, log_a, a, b)

    F = betarnd(a, b);

    if rand(1) < 0.5
        x_prop = x * F;
        log_a = log_a + log(F);
    else
        x_prop = x / F;
        log_a = log_a - log(F);
    end

end
