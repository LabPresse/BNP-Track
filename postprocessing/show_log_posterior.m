function [log_posterior, log_obs] = show_log_posterior_no_photostate(chain)
    log_posterior = zeros(length(chain.i), 1);
    log_obs = zeros(length(chain.i), 1);

    for i = 1:length(log_posterior)
        X = reshape(chain.Xm(i, :), [], chain.params.M);
        Y = reshape(chain.Ym(i, :), [], chain.params.M);
        Z = reshape(chain.Zm(i, :), [], chain.params.M);

        log_obs(i) = log_observation( ...
            chain.F(i), ...
            chain.h(i), ...
            chain.G(i), ...
            X(:, chain.bm(i, :)), ...
            Y(:, chain.bm(i, :)), ...
            Z(:, chain.bm(i, :)), ...
            chain.params ...
        );
        log_posterior(i) = log_obs(i) ...
            + log_load(chain.bm(i, :), chain.params) ...
            + log_trajectory(X, Y, Z, chain.D(i), chain.params) ...
            + log_gamma(chain.F(i), chain.params.F_prior_phi, chain.params.F_prior_psi) ...
            + log_gamma(chain.h(i), chain.params.h_prior_phi, chain.params.h_prior_psi) ...
            + log_invgamma(chain.D(i), chain.params.D_prior_phi, chain.params.D_prior_chi) ...
            + log_invgamma(chain.G(i), chain.params.G_prior_phi, chain.params.G_prior_chi);
    end

    figure;
    tiledlayout(4, 1);
    ax1 = nexttile;
    plot(chain.i(2:end), log_posterior(2:end))
    title('Log posterior')
    
    ax2 = nexttile;
    plot(chain.i(2:end), chain.D(2:end))
    title('Diffusion coefficient')
    
    ax3 = nexttile;
    hold on
    plot(chain.i(2:end), log_obs(2:end))
    hold on
    log_obs_gnd = log_observation_gnd(chain.params);
    plot([chain.i(2),chain.i(end)], [log_obs_gnd,log_obs_gnd])
    title('Log likelihood')
    
    ax4 = nexttile;
    B = sum(chain.bm, 2);
    stairs(chain.i(2:end), B(2:end))
    title('Number of active loads')

    linkaxes([ax1 ax2 ax3 ax4], 'x')

end

function p = log_observation( ...
        bkndflux, ...
        emisrate, ...
        gain, ...
        activeX, ...
        activeY, ...
        activeZ, ...
        params ...
    )
    px_PSF = zeros(params.Px, params.Py, params.N);

    for n = 1:params.N

        for m = 1:size(activeX, 2)

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

    effective_cnt = params.w_cnt / params.f / gain;
    log_eff_cnt = log(effective_cnt);

    %%
    effective_photon = params.dt_exp / params.f * (bkndflux * params.px_area + emisrate * px_PSF);
    p = sum(effective_photon .* log_eff_cnt - gammaln(effective_photon) - effective_cnt, 'all');
end

function p = log_observation_gnd(params)

    tn_bnd = params.dt_stp * (0:params.N);
    tn_min = tn_bnd(1:params.N) + 0.5 * (params.dt_stp - params.dt_exp);
    tn_max = tn_bnd(2:params.N + 1) - 0.5 * (params.dt_stp - params.dt_exp);

    dt = params.dt_exp / 1000;
    tk = (tn_bnd(1):dt:tn_bnd(end) + dt)';

    %% generate images
    u_cnt = zeros(length(params.x_bnd)-1, length(params.y_bnd)-1, params.N);

    % first gather emitter contributions
    m_ind = params.ground.Sk == 2;

    for n = 1:params.N
        k_ind = (tn_min(n) <= tk) & (tk < tn_max(n));

        for m = 1:size(params.ground.Sk,2)

            for k = find(m_ind(:, m) & k_ind)'
                u_cnt(:, :, n) = u_cnt(:, :, n) + get_g_PSF(params.ground.Xk(k, m), params.ground.Yk(k, m), params.ground.Zk(k, m), params.x_bnd, params.y_bnd, params.PSF_params);
            end

        end

    end

    u_cnt = bsxfun(@plus, params.dt_exp * params.ground.F * params.px_area, params.ground.h * dt * u_cnt);
    
    effective_cnt = params.w_cnt / params.f / params.ground.G;
    log_eff_cnt = log(effective_cnt);

    %%
    effective_photon = u_cnt / params.f;
    p = sum(effective_photon .* log_eff_cnt - gammaln(effective_photon) - effective_cnt, 'all');
end

function p = log_load(b, params)
    p = sum( ...
        -b * log1p((params.M - 1) / params.bm_prior_gamma) ...
        - (~b) * log1p(params.bm_prior_gamma / (params.M - 1)) ...
    );
end

function p = log_trajectory(X, Y, Z, D, params)
    p = -0.25 * sum(sum(diff(X).^2 + diff(Y).^2 + diff(Z).^2, 2) ./ diff(params.t_mid)) ...
        - 1.5 * (size(Z, 1) - 1) * log(D) ...
        - sum(X(1, :) - params.Xm_prior_mu).^2 / (2 * params.Xm_prior_sg^2) ...
        - sum(Y(1, :) - params.Ym_prior_mu).^2 / (2 * params.Ym_prior_sg^2) ...
        - sum(Z(1, :) - params.Zm_prior_mu).^2 / (2 * params.Zm_prior_sg^2);
end

function p = log_invgamma(x, c, d)
    p =- (c + 1) * log(x) - c * d / x;
end

function p = log_gamma(x, a, b)
    p = (a - 1) * log(x) - a * x / b;
end
