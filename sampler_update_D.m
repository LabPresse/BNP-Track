function D = sampler_update_D(active_num, activeX, activeY, activeZ, params)

    if active_num > 0
        D = params.D_prior_phi * params.D_prior_chi ...
            + 0.25 * sum(sum(diff(activeX).^2 + diff(activeY).^2 + diff(activeZ).^2, 2) ./ diff(params.t_mid));
        D = D / randg(params.D_prior_phi + 1.5 * active_num * (params.N * params.K - 1));
    else
        D = params.D_prior_phi * params.D_prior_chi;
        D = D / randg(params.D_prior_phi);
    end

end
