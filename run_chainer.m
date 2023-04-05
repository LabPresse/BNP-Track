clear
format compact

save_file = './example.mat';

if exist(save_file, 'file')
    load(save_file)
else
    %% simulate data
    [w_cnt, x_bnd, y_bnd, dt_stp, dt_exp, units, TEMP, seed] = generate_synthetic_data(50, 3);
    opts.ground = TEMP.ground;
    opts.PSF_params = TEMP.PSF_params;
    opts.f = TEMP.f;

    %% set options
    opts.x_bnd = x_bnd;
    opts.y_bnd = y_bnd;
    opts.w_cnt = w_cnt;
    opts.dt_stp = dt_stp;
    opts.dt_exp = dt_exp;
    opts.units = units;

    %% init chain
    chain = chainer_main([], 0, opts, [], []);
    chain.seed = seed;

    save(save_file, 'chain');
    disp(['SAVED: ', save_file])
end

%% expand
while true

    d_length = 10;

    chain = chainer_main(chain, d_length, [], true, true);

    if chain.sizeGB > 0.5
        chain = chainer_main(chain, -fix(chain.length / 2), [], true, true);
    end

    save(save_file)
    disp(['SAVED: ', save_file])

    %figure(99)
    %show_emitter_tracks(chain)

end
