function w_cnt = generate_synthetic_measurments(u_cnt, gain, f, ...
        units, win, ...
        x_bnd, y_bnd)

    %% add shot noise
    w_cnt = f * gain * randg(u_cnt / f);

    %%
    if ~isempty(win)

        figure(win)
        set(gcf, 'name', 'Measurments')
        clf
        colormap('gray')

        N = size(w_cnt, 3);

        u_lim = [0 max(u_cnt(:))];
        w_lim = [0 max(w_cnt(:))];

        n = 1;

        subplot(6, 3, [1 7])
        plot_tiles(x_bnd, y_bnd, u_cnt(:, :, n))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        ylabel(['y (', units.length, ')'])
        title(['n=', num2str(n), ' - N=', num2str(N)])
        caxis(u_lim);

        subplot(6, 3, [10 16])
        plot_tiles(x_bnd, y_bnd, w_cnt(:, :, n))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        ylabel(['y (', units.length, ')'])
        xlabel(['x (', units.length, ')'])
        caxis(w_lim);

        n = N;

        subplot(6, 3, [2 8])
        plot_tiles(x_bnd, y_bnd, u_cnt(:, :, n))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        title(['n=', num2str(n), ' - N=', num2str(N)])
        caxis(u_lim);

        subplot(6, 3, [11 17])
        plot_tiles(x_bnd, y_bnd, w_cnt(:, :, n))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        xlabel(['x (', units.length, ')'])
        caxis(w_lim);

        subplot(6, 3, [3 9])
        plot_tiles(x_bnd, y_bnd, mean(u_cnt, 3))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        title('mean')
        caxis(u_lim);

        subplot(6, 3, [12 18])
        plot_tiles(x_bnd, y_bnd, mean(w_cnt, 3))
        axis image
        xlim([x_bnd(1) x_bnd(end)])
        ylim([y_bnd(1) y_bnd(end)])
        xlabel(['x (', units.length, ')'])
        caxis(w_lim);

    end

    %%
    if nargout < 1
        clear w_cnt
    end
