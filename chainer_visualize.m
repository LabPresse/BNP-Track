function Gim = chainer_visualize(Gim, chain)

    % mask_sm = nan(chain.params.N,chain.params.M);
    % mask_sm(chain.sample.sm==2) = 1;

    mask_Sm = get_mask_Sm(chain.sample.sm, chain.params);

    Bm = sum(chain.bm, 2);
    Bm(Bm == 0) = nan;

    %% init
    if isempty(Gim)

        i = double(chain.i(1)) + chain.params.i_skip * (0:chain.length - 1)';

        mum = 3;
        num = 7;

        figure(10)
        set(gcf, 'windowstyle', 'docked')
        clf

        %----------------------------------------------------------------------
        subplot(num, mum, 0 * mum + [1 mum - 1]);
        axis off
        title('Current sample')

        %----------------------------------------------------------------------
        subplot(num, mum, 1 * mum + [1 2 * mum - 1]);
        Gim.ax_SAMP_x = plot(chain.params.t_mid, chain.sample.Xm, 'b', ...
            chain.params.t_mid, chain.sample.Xm .* mask_Sm, '.-');

        if isfield(chain.params, 'ground')
            temp = line( ...
                chain.params.ground.tk, ...
                chain.params.ground.Xk, ...
                'color', 'r', 'LineWidth', 2 ...
            );

            for m = 1:length(temp)
                uistack(temp(m), 'bottom');
            end

        end

        clear temp
        xlim(chain.params.t_bnd([1 end]))

        for m = 1:chain.params.M

            if chain.sample.bm(m)
                Gim.ax_SAMP_x(m).LineStyle = '-';
                Gim.ax_SAMP_x(m).LineWidth = 2;
                Gim.ax_SAMP_x(m + chain.params.M).Color = 'g';
            else
                Gim.ax_SAMP_x(m).LineStyle = ':';
                Gim.ax_SAMP_x(m + chain.params.M).Color = 'b';
            end

        end

        add_prior(gca, 'norm', chain.params.Xm_prior_mu, chain.params.Xm_prior_sg)
        ylabel(['X^m (', chain.params.units.length, ')'])
        box off
        grid on
        line(chain.params.t_bnd(end - 2) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_bnd(end - 1) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_mid(chain.params.t_idx(end - 1, :), :) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', ':');
        line(get(gca, 'xlim'), chain.params.x_bnd(1) * [1 1], 'color', 'k', 'linestyle', '--');
        line(get(gca, 'xlim'), chain.params.x_bnd(end) * [1 1], 'color', 'k', 'linestyle', '--');
        ylim(get(gca, 'Ylim'))

        %----------------------------------------------------------------------
        subplot(num, mum, 3 * mum + [1 2 * mum - 1]);
        Gim.ax_SAMP_y = plot(chain.params.t_mid, chain.sample.Ym, 'b', ...
            chain.params.t_mid, chain.sample.Ym .* mask_Sm, '.-');

        if isfield(chain.params, 'ground')
            temp = line( ...
                chain.params.ground.tk, ...
                chain.params.ground.Yk, ...
                'color', 'r', 'LineWidth', 2 ...
            );

            for m = 1:length(temp)
                uistack(temp(m), 'bottom');
            end

        end

        clear temp
        xlim(chain.params.t_bnd([1 end]))

        for m = 1:chain.params.M

            if chain.sample.bm(m)
                Gim.ax_SAMP_y(m).LineStyle = '-';
                Gim.ax_SAMP_y(m).LineWidth = 2;
                Gim.ax_SAMP_y(m + chain.params.M).Color = 'g';
            else
                Gim.ax_SAMP_y(m).LineStyle = ':';
                Gim.ax_SAMP_y(m + chain.params.M).Color = 'b';
            end

        end

        add_prior(gca, 'norm', chain.params.Ym_prior_mu, chain.params.Ym_prior_sg)
        ylabel(['Y^m (', chain.params.units.length, ')'])
        box off
        grid on
        line(chain.params.t_bnd(end - 2) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_bnd(end - 1) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_mid(chain.params.t_idx(end - 1, :), :) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', ':');
        line(get(gca, 'xlim'), chain.params.y_bnd(1) * [1 1], 'color', 'k', 'linestyle', '--');
        line(get(gca, 'xlim'), chain.params.y_bnd(end) * [1 1], 'color', 'k', 'linestyle', '--');
        ylim(get(gca, 'Ylim'))

        %----------------------------------------------------------------------
        subplot(num, mum, 5 * mum + [1 2 * mum - 1]);
        Gim.ax_SAMP_z = plot(chain.params.t_mid, chain.sample.Zm, 'b', ...
            chain.params.t_mid, chain.sample.Zm .* mask_Sm, '.-');

        if isfield(chain.params, 'ground')
            temp = line( ...
                chain.params.ground.tk, ...
                chain.params.ground.Zk, ...
                'color', 'r', 'LineWidth', 2 ...
            );

            for m = 1:length(temp)
                uistack(temp(m), 'bottom');
            end

        end

        clear temp
        xlim(chain.params.t_bnd([1 end]))

        for m = 1:chain.params.M

            if chain.sample.bm(m)
                Gim.ax_SAMP_z(m).LineStyle = '-';
                Gim.ax_SAMP_z(m).LineWidth = 2;
                Gim.ax_SAMP_z(m + chain.params.M).Color = 'g';
            else
                Gim.ax_SAMP_z(m).LineStyle = ':';
                Gim.ax_SAMP_z(m + chain.params.M).Color = 'b';
            end

        end

        add_prior(gca, 'symnorm', chain.params.Zm_prior_mu, chain.params.Zm_prior_sg)
        ylabel(['z^m (', chain.params.units.length, ')'])
        box off
        grid on
        line(chain.params.t_bnd(end - 2) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_bnd(end - 1) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', '-');
        line(chain.params.t_mid(chain.params.t_idx(end - 1, :), :) * [1 1], get(gca, 'ylim'), 'color', 'k', 'linestyle', ':');
        %     line(get(gca,'xlim'),chain.params.x_bnd(  1)*[1 1],'color','k','linestyle','--');
        %     line(get(gca,'xlim'),chain.params.x_bnd(end)*[1 1],'color','k','linestyle','--');
        line(get(gca, 'xlim'), 0 * [1 1], 'color', 'k', 'linestyle', '--');
        line(get(gca, 'xlim'), +chain.params.PSF_params.z_ref * [1 1], 'color', 'k', 'linestyle', ':');
        line(get(gca, 'xlim'), -chain.params.PSF_params.z_ref * [1 1], 'color', 'k', 'linestyle', ':');
        ylim(get(gca, 'Ylim'))

        % add PSF
        z_temp = linspace(min(get(gca, 'Ylim')), max(get(gca, 'Ylim')));
        p_temp = min(get(gca, 'Xlim')) + 0.15 * (max(get(gca, 'Xlim')) - min(get(gca, 'Xlim'))) ./ (1 + (z_temp / chain.params.PSF_params.z_ref) .^ 2);
        line(p_temp, z_temp, 'color', 'k', 'linestyle', ':');

        xlabel(['t (', chain.params.units.time, ')'])

        %----------------------------------------------------------------------
        subplot(num, mum, [1 * mum 1 * mum]);
        xlim([max(0, i(1) - 1), i(end) + 1])
        ylim([0 min(chain.params.M, 2 * nnz(chain.sample.bm)) + 1])
        set(gca, 'YTick', 1:chain.params.M)

        for m = 1:chain.params.M
            Gim.ax_MCMC_b{1}(m) = line(i, get_a(chain.bm(:, m), m), 'marker', '.');
        end

        Gim.ax_MCMC_b{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        Gim.ax_MCMC_b{3} = line(i, Bm, 'color', sum([0 0; 0 1] * lines(2)));

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), chain.params.ground.M * [1 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        ylabel('b^m')
        title('MCMC')
        set(gca, 'Yaxislocation', 'right')
        box off
        grid on

        %----------------------------------------------------------------------
        subplot(num, mum, [2 * mum 3 * mum]);
        Gim.ax_MCMC_D{1} = plot(i, chain.D, '.-');
        xlim([max(0, i(1) - 1), i(end) + 1])

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), chain.params.ground.D * [1 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        % ylim([0 chain.params.M+1])
        Gim.ax_MCMC_D{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        ylabel({'D'; ['(', chain.params.units.length, '^2/', chain.params.units.time, ')']})
        set(gca, 'Yaxislocation', 'right')
        add_prior(gca, 'invgamma', chain.params.D_prior_phi, chain.params.D_prior_chi)
        box off
        grid on

        %----------------------------------------------------------------------
        subplot(num, mum, [4 * mum 4 * mum]);
        Gim.ax_MCMC_P{1} = plot(i, ones(length(chain.D), chain.params.M), '.-');
        %     set(gca,'YScale','Log')
        uistack(Gim.ax_MCMC_P{1}(3), 'bottom')
        xlim([max(0, i(1) - 1), i(end) + 1])

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), [1; 1] * [1, 1, 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        % ylim([0 chain.params.M+1])
        Gim.ax_MCMC_P{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        ylabel({'\pi_{\sigma}'; '(1/stp)'})
        set(gca, 'Yaxislocation', 'right')
        %add_prior(gca, 'beta', chain.params.Pj_prior_A(1), chain.params.Pj_prior_B(1))
        %add_prior(gca, 'beta', chain.params.Pj_prior_A(2), chain.params.Pj_prior_B(2))
        %add_prior(gca, 'beta', chain.params.Pj_prior_A(3), chain.params.Pj_prior_B(3))
        box off
        grid on

        %----------------------------------------------------------------------
        subplot(num, mum, [5 * mum 5 * mum]);
        Gim.ax_MCMC_F{1} = plot(i, chain.F, '.-');
        xlim([max(0, i(1) - 1), i(end) + 1])

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), chain.params.ground.F * [1 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        % ylim([0 chain.params.M+1])
        Gim.ax_MCMC_F{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        ylabel({'F'; ['(\gamma/', chain.params.units.time, '/', chain.params.units.length, '^2)']})
        set(gca, 'Yaxislocation', 'right')
        add_prior(gca, 'gamma', chain.params.F_prior_phi, chain.params.F_prior_psi)
        box off
        grid on

        %----------------------------------------------------------------------
        subplot(num, mum, [6 * mum 6 * mum]);
        Gim.ax_MCMC_h{1} = plot(i, chain.h, '.-');
        xlim([max(0, i(1) - 1), i(end) + 1])

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), chain.params.ground.h * [1 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        % ylim([0 chain.params.M+1])
        Gim.ax_MCMC_h{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        ylabel({'h'; ['(\gamma/', chain.params.units.time, ')']})
        set(gca, 'Yaxislocation', 'right')
        add_prior(gca, 'gamma', chain.params.h_prior_phi, chain.params.h_prior_psi)
        box off
        grid on

        %----------------------------------------------------------------------
        subplot(num, mum, [7 * mum 7 * mum]);
        Gim.ax_MCMC_G{1} = plot(i, chain.G, '.-');
        xlim([max(0, i(1) - 1), i(end) + 1])

        if isfield(chain.params, 'ground')
            temp = line(get(gca, 'XLim'), chain.params.ground.G * [1 1], 'color', 'r', 'linewidth', 2);
            uistack(temp, 'top');
        end

        clear temp
        % ylim([0 chain.params.M+1])
        Gim.ax_MCMC_G{2} = line(chain.sample.i * [1 1], get(gca, 'Ylim'), 'color', 'k', 'linestyle', ':');
        ylabel({'G'; ['(', chain.params.units.image, '/\gamma)']})
        set(gca, 'Yaxislocation', 'right')
        add_prior(gca, 'invgamma', chain.params.G_prior_phi, chain.params.G_prior_chi)
        box off
        grid on

        xlabel('MCMC iteration (i)')

    end

    for m = 1:chain.params.M
        Gim.ax_SAMP_x(m).YData = chain.sample.Xm(:, m);
        Gim.ax_SAMP_x(m + chain.params.M).YData = chain.sample.Xm(:, m) .* mask_Sm(:, m);

        Gim.ax_SAMP_y(m).YData = chain.sample.Ym(:, m);
        Gim.ax_SAMP_y(m + chain.params.M).YData = chain.sample.Ym(:, m) .* mask_Sm(:, m);

        Gim.ax_SAMP_z(m).YData = chain.sample.Zm(:, m);
        Gim.ax_SAMP_z(m + chain.params.M).YData = chain.sample.Zm(:, m) .* mask_Sm(:, m);

        if chain.sample.bm(m)
            Gim.ax_SAMP_x(m).LineStyle = '-';
            Gim.ax_SAMP_x(m + chain.params.M).Color = 'g';

            Gim.ax_SAMP_y(m).LineStyle = '-';
            Gim.ax_SAMP_y(m + chain.params.M).Color = 'g';

            Gim.ax_SAMP_z(m).LineStyle = '-';
            Gim.ax_SAMP_z(m + chain.params.M).Color = 'g';
        else
            Gim.ax_SAMP_x(m).LineStyle = ':';
            Gim.ax_SAMP_x(m + chain.params.M).Color = 'b';

            Gim.ax_SAMP_y(m).LineStyle = ':';
            Gim.ax_SAMP_y(m + chain.params.M).Color = 'b';

            Gim.ax_SAMP_z(m).LineStyle = ':';
            Gim.ax_SAMP_z(m + chain.params.M).Color = 'b';
        end

        Gim.ax_MCMC_b{1}(m).YData = get_a(chain.bm(:, m), m);
    end

    Gim.ax_MCMC_D{1}.YData = chain.D;
    Gim.ax_MCMC_F{1}.YData = chain.F;
    Gim.ax_MCMC_h{1}.YData = chain.h;
    Gim.ax_MCMC_G{1}.YData = chain.G;
    Gim.ax_MCMC_P{1}(1).YData = ones(length(chain.D), 1);
    Gim.ax_MCMC_P{1}(2).YData = ones(length(chain.D), 1);
    Gim.ax_MCMC_P{1}(3).YData = ones(length(chain.D), 1);

    Gim.ax_MCMC_b{2}.XData = chain.sample.i * [1 1];
    Gim.ax_MCMC_D{2}.XData = chain.sample.i * [1 1];
    Gim.ax_MCMC_F{2}.XData = chain.sample.i * [1 1];
    Gim.ax_MCMC_h{2}.XData = chain.sample.i * [1 1];
    Gim.ax_MCMC_G{2}.XData = chain.sample.i * [1 1];
    Gim.ax_MCMC_P{2}.XData = chain.sample.i * [1 1];

    Gim.ax_MCMC_b{3}.YData = Bm;

    drawnow
end

%%
function a = get_a(b, m)

    a = nan(length(b), 1);
    a(b) = m;
end