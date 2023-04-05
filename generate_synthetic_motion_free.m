function [xmk, ymk, zmk] = generate_synthetic_motion_free(dt, K, M, D, ...
        x0_mu, x0_sg, ...
        y0_mu, y0_sg, ...
        z0_mu, z0_sg, ...
        units, win, tk)

    %% generate initial positions
    xm0 = x0_sg * randn(1, M) + x0_mu;
    ym0 = y0_sg * randn(1, M) + y0_mu;
    zm0 = z0_sg * randn(1, M) + z0_mu * (2 * (double(rand(1, M) < 0.5) - 0.5));

    %% generate free Brownian motion
    s = sqrt(2 * D * dt);

    xmk = cumsum([xm0; s * randn(K - 1, M)]);
    ymk = cumsum([ym0; s * randn(K - 1, M)]);
    zmk = cumsum([zm0; s * randn(K - 1, M)]);

    %% plot the simulated trajectories
    figure(win)
    set(gcf, 'name', 'Motion')
    clf

    subplot(3, 1, 1)
    plot(tk, xmk, 'g');
    xlim([0 tk(end) + 0.5 * dt])
    ylabel(['x^m_k (', units.length, ')'])
    box off

    subplot(3, 1, 2)
    plot(tk, ymk, 'g');
    xlim([0 tk(end) + 0.5 * dt])
    ylabel(['y^m_k (', units.length, ')'])
    box off

    subplot(3, 1, 3)
    plot(tk, zmk, 'g');
    xlim([0 tk(end) + 0.5 * dt])
    ylabel(['z^m_k (', units.length, ')'])
    box off

    xlabel(['t_k (', units.time, ')'])

end
