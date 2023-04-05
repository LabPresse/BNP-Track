function u_cnt = generate_synthetic_images(N, M, x_bnd, y_bnd, ...
        dt_exp, ...
        PSF_params, F, h, ...
        dt, tk, xmk, ymk, zmk, ...
        units, win, ...
        dt_stp, tn_bnd)

    %% get sizes
    Px = length(x_bnd) - 1; % num of pixels in x-coord
    Py = length(y_bnd) - 1; % num of pixels in y-coord

    %% ensure correct shapes
    x_bnd = reshape(x_bnd, Px + 1, 1); % coloumn
    y_bnd = reshape(y_bnd, 1, Py + 1); % row

    %% establish bounds for exposures
    exp_bnd = tn_bnd(1:end - 1) + 0.5 * (dt_stp - dt_exp);
    exp_bnd(2, :) = exp_bnd(1, :) + dt_exp;

    %% generate images
    u_cnt = zeros(Px, Py, N);

    % first gather emitter contributions
    for n = 1:N
        k_ind = (exp_bnd(1, n) <= tk) & (tk < exp_bnd(2, n));

        for m = 1:M

            for k = find(k_ind)'
                u_cnt(:, :, n) = u_cnt(:, :, n) + get_g_PSF(xmk(k, m), ymk(k, m), zmk(k, m), x_bnd, y_bnd, PSF_params);
            end

        end

    end

    u_cnt = dt_exp * F * diff(x_bnd) * diff(y_bnd) + h * dt * u_cnt;

    %% demo
    figure(win)
    set(gcf, 'name', 'Images')
    clf
    colormap('gray')

    n = 1;

    k_ind = (tn_bnd(n) <= tk) & (tk < tn_bnd(n + 1));
    m_ind = (exp_bnd(1, n) <= tk) & (tk < exp_bnd(2, n));

    subplot(6, 3, [1 7])
    plot_tiles(x_bnd, y_bnd, u_cnt(:, :, n))
    axis image
    xlim([x_bnd(1) x_bnd(end)])
    ylim([y_bnd(1) y_bnd(end)])
    ylabel(['y (', units.length, ')'])
    title(['n=', num2str(n), ' - N=', num2str(N)])

    subplot(6, 3, [10 16])
    plot_tiles(x_bnd, y_bnd, u_cnt(:, :, n))
    axis image
    xlim([x_bnd(1) x_bnd(end)])
    ylim([y_bnd(1) y_bnd(end)])

    for m = 1:M
        line(xmk(k_ind, m), ymk(k_ind, m), 'linestyle', '-', 'color', 'b');
    end

    for m = 1:M
        line(xmk(m_ind, m), ymk(m_ind, m), 'linestyle', '-', 'color', 'g');
    end

    xlabel(['x (', units.length, ')'])
    ylabel(['y (', units.length, ')'])

    subplot(6, 3, [2 6])
    plot(tk(k_ind), xmk(k_ind, :), 'b');

    for m = 1:M
        line(tk(m_ind), xmk(m_ind, m), 'linestyle', '-', 'color', 'g', 'marker', '.');
    end

    xlim(0.5 * (exp_bnd(1, n) + exp_bnd(2, n)) + 0.6 * dt_stp * [-1 +1])
    line(tn_bnd(n) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(tn_bnd(n + 1) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(exp_bnd(1, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    line(exp_bnd(2, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    ylabel(['x^m_k (', units.length, ')'])
    box off
    p1 = patch([tn_bnd(n) exp_bnd(1, n) exp_bnd(1, n) tn_bnd(n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p2 = patch([tn_bnd(n + 1) exp_bnd(2, n) exp_bnd(2, n) tn_bnd(n + 1)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p3 = patch([exp_bnd(1, n) exp_bnd(2, n) exp_bnd(2, n) exp_bnd(1, n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.75 0.75 0.75], 'edgecolor', [0.75 0.75 0.75]);
    uistack(p3, 'bottom')
    uistack(p1, 'bottom')
    uistack(p2, 'bottom')

    title('Exposure during one aquisition period')

    subplot(6, 3, [8 12])
    plot(tk(k_ind), ymk(k_ind, :), 'b');

    for m = 1:M
        line(tk(m_ind), ymk(m_ind, m), 'linestyle', '-', 'color', 'g', 'marker', '.');
    end

    xlim(0.5 * (exp_bnd(1, n) + exp_bnd(2, n)) + 0.6 * dt_stp * [-1 +1])
    line(tn_bnd(n) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(tn_bnd(n + 1) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(exp_bnd(1, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    line(exp_bnd(2, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    ylabel(['y^m_k (', units.length, ')'])
    box off
    p1 = patch([tn_bnd(n) exp_bnd(1, n) exp_bnd(1, n) tn_bnd(n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p2 = patch([tn_bnd(n + 1) exp_bnd(2, n) exp_bnd(2, n) tn_bnd(n + 1)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p3 = patch([exp_bnd(1, n) exp_bnd(2, n) exp_bnd(2, n) exp_bnd(1, n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.75 0.75 0.75], 'edgecolor', [0.75 0.75 0.75]);
    uistack(p3, 'bottom')
    uistack(p1, 'bottom')
    uistack(p2, 'bottom')

    subplot(6, 3, [14 18])
    plot(tk(k_ind), zmk(k_ind, :), 'b');

    for m = 1:M
        line(tk(m_ind), zmk(m_ind, m), 'linestyle', '-', 'color', 'g', 'marker', '.');
    end

    xlim(0.5 * (exp_bnd(1, n) + exp_bnd(2, n)) + 0.6 * dt_stp * [-1 +1])
    line(tn_bnd(n) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(tn_bnd(n + 1) * [1 1], get(gca, 'ylim'), 'color', 'k')
    line(exp_bnd(1, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    line(exp_bnd(2, n) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '--')
    line(get(gca, 'Xlim'), 0 * [1 1], 'color', 'k', 'linestyle', ':')
    ylabel(['z^m_k (', units.length, ')'])
    xlabel(['t_k (', units.time, ')'])
    box off
    p1 = patch([tn_bnd(n) exp_bnd(1, n) exp_bnd(1, n) tn_bnd(n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p2 = patch([tn_bnd(n + 1) exp_bnd(2, n) exp_bnd(2, n) tn_bnd(n + 1)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.50 0.50 0.50], 'edgecolor', [0.50 0.50 0.50]);
    p3 = patch([exp_bnd(1, n) exp_bnd(2, n) exp_bnd(2, n) exp_bnd(1, n)], [1 1] * diag(get(gca, 'ylim')) * [1 1 0 0; 0 0 1 1], [0.75 0.75 0.75], 'edgecolor', [0.75 0.75 0.75]);
    uistack(p3, 'bottom')
    uistack(p1, 'bottom')
    uistack(p2, 'bottom')

end
