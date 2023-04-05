function [emisrate, loads, rec] = ...
        sampler_update_rates_and_loads( ...
        emisrate, ...
        loads, ...
        bkndflux, ...
        x, ...
        y, ...
        z, ...
        gain, ...
        rec, ...
        T, ...
        params, ...
        batch_size ...
    )

    emisrate_prop = get_proposal(emisrate, params.MH_sc(2), 1);
    % emisrate_prop = params.ground.h;

    emitter_contr = zeros(params.P, params.N, params.M);

    bknd_contr = params.dt_exp / params.f * bkndflux * params.px_area_col;
    log_signal = params.w_cnt_log - log(gain);
    log_prior = [log(params.M - 1); log(params.bm_prior_gamma)];

    active_idx = find(loads);
    inactive_idx = setdiff(1:params.M, active_idx);
    shuffled_idx1 = active_idx(randperm(length(active_idx)));
    shuffled_idx2 = inactive_idx(randperm(length(inactive_idx)));
    sample_order = [shuffled_idx1, shuffled_idx2];

    batch = sample_order(1:batch_size);
    active_and_out = setdiff(find(loads), batch);

    for m = 1:params.M

        for k = 1:params.K
            % params.ak saves the unitless coefficients of the composite trapezoidal quadrature
            emitter_contr(:, :, m) = emitter_contr(:, :, m) + reshape( ...
            params.ak(k) * get_g_PSF_2( ...
                x(params.t_idx(:, k), m), ...
                y(params.t_idx(:, k), m), ...
                z(params.t_idx(:, k), m), ...
                params.x_bnd, ...
                params.y_bnd, ...
                params.PSF_params), params.P, params.N, 1);
        end

    end

    emitter_contr = params.dt_exp / params.f * emisrate_prop * emitter_contr;

    base = bknd_contr + sum(emitter_contr(:, :, active_and_out), 3);

    log_post_prop = get_log_post(emitter_contr(:, :, batch), base, log_signal, log_prior, T);
    emitter_contr = emitter_contr / emisrate_prop * emisrate;
    log_post = get_log_post(emitter_contr(:, :, batch), base, log_signal, log_prior, T);

    log_r = log_sum_exp(log_post_prop) - log_sum_exp(log_post) ...
        + (params.h_prior_phi -1) * log(emisrate_prop / emisrate) ...
        + params.h_prior_phi * (emisrate - emisrate_prop) / params.h_prior_psi;

    if log(rand) < log_r
        emisrate = emisrate_prop;
        log_post = log_post_prop;
        emitter_contr = emitter_contr * emisrate_prop / emisrate;
        rec(1) = rec(1) + 1;
    end

    rec(2) = rec(2) + 1;

    [~, result] = max(log_post - log(-log(rand(length(log_post), 1))));
    loads(batch) = logical(dec2bin(result - 1, batch_size) - 48);

    for idx = batch_size + 1:batch_size:params.M
        batch = sample_order(idx:min(idx + batch_size - 1, params.M));
        active_and_out = setdiff(find(loads), batch);
        base = bknd_contr + sum(emitter_contr(:, :, active_and_out), 3);
        log_post = get_log_post(emitter_contr(:, :, batch), base, log_signal, log_prior, T);
        [~, result] = max(log_post - log(-log(rand(length(log_post), 1))));
        loads(batch) = logical(dec2bin(result - 1, batch_size) - 48);
    end

end

function x_prop = get_proposal(x, a, b)

    F = betarnd(a, b);

    if rand(1) < 0.5
        x_prop = x * F;
    else
        x_prop = x / F;
    end

end

function log_post = get_log_post(batch_contr, base, log_signal, log_prior, T)
    batch_size = size(batch_contr, 3);
    log_post = zeros(2^batch_size, 1);

    for m = 1:length(log_post)
        new_on = logical(dec2bin(m - 1, batch_size) - 48);
        emitter_contr = base + sum(batch_contr(:, :, new_on), 3);
        log_post(m) = sum(emitter_contr .* log_signal - gammaln(emitter_contr), 'all') / T...
            + nnz(new_on) * log_prior(2) + (batch_size - nnz(new_on)) * log_prior(1);
    end

end

function log_x = log_sum_exp(log_x)

    log_x = [-inf, sort(reshape(log_x(log_x > -inf), 1, []))];

    while length(log_x) > 1
        log_x = [max(log_x(1), log_x(2)) + log1p(exp(-abs(log_x(1) - log_x(2)))), log_x(3:end)];
    end

end
