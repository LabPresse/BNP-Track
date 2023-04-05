function mask_Sm = get_mask_Sm(sm, params)

    mask_Sm = nan(params.N * params.K, params.M);

    for m = 1:params.M
        mask_Sm(params.t_idx(sm(:, m) == 2, :), m) = 1;
    end

    %mask_sm(chain.sample.sm==2) = 1;
