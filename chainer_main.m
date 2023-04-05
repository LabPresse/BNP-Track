function chain = chainer_main(chain_init, d_length, opts, flag_status, flag_visual)
    % to init:
    % chain = chainer_main([], 0, opts, true, []);
    % to expand:
    % chain = chainer_main(chain, +25, [], true, true);
    % to reduce:
    % chain = chainer_main(chain, -10, [], true, []);

    rng('shuffle')

    % init chain --------------------------------------------------------------
    if d_length == 0
        tic_id = tic;

        % MCMC
        chain.params = chainer_init_params(opts);
        chain.length = 1;
        chain.ledger = nan(0, 2);
        chain.sizeGB = nan;
        chain.record = [];
        chain.sample = [];

        chain.sample = chainer_init_sample(chain.params);

        % history
        chain.i = cast(chain.sample.i, 'uint64');
        chain.T = cast(chain.sample.T, 'double');
        chain.bm = cast(chain.sample.bm, 'logical');
        chain.Xm = cast(chain.sample.Xm(:)', 'single'); % reshape
        chain.Ym = cast(chain.sample.Ym(:)', 'single'); % reshape
        chain.Zm = cast(chain.sample.Zm(:)', 'single'); % reshape
        chain.D = cast(chain.sample.D, 'double');
        chain.F = cast(chain.sample.F, 'double');
        chain.h = cast(chain.sample.h, 'double');
        chain.G = cast(chain.sample.G, 'double');

        if flag_status
            disp(['CHAINER: chain initiated (total time = ', num2str(toc(tic_id)), ' s)'])
        end

        % expand chain ------------------------------------------------------------
    elseif d_length > 0
        tic_id = tic;

        chain.params = chain_init.params;
        chain.length = chain_init.length + d_length;
        chain.ledger = chain_init.ledger;
        chain.sizeGB = nan;
        chain.record = chain_init.record;
        chain.sample = chain_init.sample;

        L = chain.params.M * chain.params.N;

        chain.i = [chain_init.i; zeros(d_length, 1, 'like', chain_init.i)];
        chain.T = [chain_init.T; zeros(d_length, 1, 'like', chain_init.T)];
        chain.bm = [chain_init.bm; false(d_length, chain.params.M)];
        chain.Xm = [chain_init.Xm; nan(d_length, L * chain.params.K, 'like', chain_init.Xm)];
        chain.Ym = [chain_init.Ym; nan(d_length, L * chain.params.K, 'like', chain_init.Ym)];
        chain.Zm = [chain_init.Zm; nan(d_length, L * chain.params.K, 'like', chain_init.Zm)];
        chain.D = [chain_init.D; nan(d_length, 1, 'like', chain_init.D)];
        chain.F = [chain_init.F; nan(d_length, 1, 'like', chain_init.F)];
        chain.h = [chain_init.h; nan(d_length, 1, 'like', chain_init.h)];
        chain.G = [chain_init.G; nan(d_length, 1, 'like', chain_init.G)];

        if flag_visual
            Gim = chainer_visualize([], chain);
        end

        %---------------------------- expand chain
        r = chain_init.length + 1;

        while r <= chain.length

            [chain.sample, chain.params] = sampler_update(chain.sample, chain.params);

            [ ...
                 chain.params, ...
                 chain.record, ...
                 chain.sample.rec_Fh, ...
                 chain.sample.rec_XYZ ...
             ] = sampler_adapt_proposals( ...
                chain.params, ...
                chain.record, ...
                chain.sample.rec_Fh, ...
                chain.sample.rec_XYZ, ...
                chain.sample.i, ...
                chain.T, ...
                chain.i ...
            );

            if mod(chain.sample.i, chain.params.i_skip) == 0

                chain.i(r) = chain.sample.i;
                chain.bm(r, :) = chain.sample.bm;
                chain.Xm(r, :) = chain.sample.Xm(:);
                chain.Ym(r, :) = chain.sample.Ym(:);
                chain.Zm(r, :) = chain.sample.Zm(:);
                chain.D(r) = chain.sample.D;
                chain.F(r) = chain.sample.F;
                chain.h(r) = chain.sample.h;
                chain.G(r) = chain.sample.G;
                chain.T(r) = chain.sample.T;

                if flag_visual
                    chainer_visualize(Gim, chain);
                end

                if flag_status
                    disp([ ...
                              'i = ', num2str(chain.sample.i, '%d'), ...
                              ' - T = ', num2str(chain.sample.T, '%#6.2f'), ...
                              ' - B = ', num2str(sum(chain.sample.bm), '%d'), ...
                              ' - D = ', num2str(chain.sample.D, '%#6.2e'), ...
                              ' - F = ', num2str(chain.sample.F, '%#6.2e'), ...
                              ' - G = ', num2str(chain.sample.G, '%#6.2e'), ...
                              ' - h = ', num2str(chain.sample.h, '%#6.2e'), ...
                              ' - rate_accep = ', num2str(chain.sample.rec_Fh(1, :) ./ chain.sample.rec_Fh(2, :) * 100, '%#6.2f'), ...
                              ' - traj_accep = ', num2str(chain.sample.rec_XYZ(1, :) ./ chain.sample.rec_XYZ(2, :) * 100, '%#6.2f'), ' %', ...
                          ])
                end

                r = r + 1;
            end

        end

        % expansion ledger
        wall_time = toc(tic_id);
        chain.ledger = [chain.ledger; double(chain.i(end)), wall_time];

        if flag_status
            disp([ ...
                      'CHAINER: chain expanded ( wall time = ', num2str(wall_time), ...
                      ' s, overall = ', num2str(sum(chain.ledger(:, 2)) / chain.ledger(end, 1)), ' s / i)' ...
                  ])
        end

        % reduce chain ------------------------------------------------------------
    elseif d_length < 0

        d_length = -d_length;

        chain.params = chain_init.params;
        chain.length = d_length;
        chain.ledger = chain_init.ledger;
        chain.sizeGB = nan;
        chain.record = chain_init.record;
        chain.sample = chain_init.sample;

        ind = mod(chain_init.length, d_length) ...
            + (floor(chain_init.length / d_length) * (1:d_length));

        chain.i = chain_init.i(ind);
        chain.T = chain_init.T(ind);
        chain.bm = chain_init.bm(ind, :);
        chain.Xm = chain_init.Xm(ind, :);
        chain.Ym = chain_init.Ym(ind, :);
        chain.Zm = chain_init.Zm(ind, :);
        chain.D = chain_init.D(ind);
        chain.F = chain_init.F(ind);
        chain.h = chain_init.h(ind);
        chain.G = chain_init.G(ind);

        chain.params.i_skip = double(chain.i(2) - chain.i(1));

        if flag_status
            disp('CHAINER: chain reduced')
        end

    end

    %% book-keeping
    chain.sizeGB = get_sizeGB(chain); % mem size

end

%% auxiliary functions

function sizeGB = get_sizeGB(chain)
    sizeGB = whos(inputname(1));
    sizeGB = sizeGB.bytes / 1024 ^ 3;
end

function str = get_int_type(L)

    if L <= intmax('uint8')
        str = 'uint8';
    elseif L <= intmax('uint16')
        str = 'uint16';
    elseif L <= intmax('uint32')
        str = 'uint32';
    elseif L <= intmax('uint64')
        str = 'uint64';
    else
        error('L way too large !!!');
    end

end

%% adapter
function [params, record, rec_Fh, rec_XYZ] = ...
        sampler_adapt_proposals(params, record, rec_Fh, rec_XYZ, i, T_chain, i_chain)

    T_idx = find(T_chain <= 2 & T_chain ~= 0, 1);

    if isempty(T_idx)
        return
    else
        i_T = i_chain(T_idx);
    end

    acc_target_Fh = [50 50 35] / 100;
    acc_target_XYZ = [25 25] / 100;
    rep = 25;

    if isempty(record)
        record = [rec_Fh, rec_XYZ];
    else
        record = record + [rec_Fh, rec_XYZ];
    end

    if i - i_T <= rep * 20

        if mod(i, rep) == 0
            disp('=== ADAPTER ===')

            acc_batch_Fh = rec_Fh(1, 1:2) ./ rec_Fh(2, 1:2);
            params.MH_sc(1) = params.MH_sc(1) ...
                / min(max(acc_batch_Fh(1) / acc_target_Fh(1), 0.1), 10);
            params.MH_sc(2) = params.MH_sc(2) ...
                / min(max(acc_batch_Fh(2) / acc_target_Fh(2), 0.1), 10);

            acc_batch_XYZ = rec_XYZ(1, 1:2) ./ rec_XYZ(2, 1:2);
            params.MH_sc(4) = params.MH_sc(4) ...
                * min(max(acc_batch_XYZ(1) / acc_target_XYZ(1), 0.1), 10);
            params.MH_sc(5) = params.MH_sc(5) ...
                * min(max(acc_batch_XYZ(2) / acc_target_XYZ(2), 0.1), 10);

            rec_Fh(:) = 0;
            rec_XYZ(:) = 0;
        end

    end

end
