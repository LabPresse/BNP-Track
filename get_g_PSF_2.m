function g = get_g_PSF_weiqing(x, y, z, x_bnd, y_bnd, PSF_params)
    % g = int_Ap G(x,y,z;x,y) dx dt
    % size(g) = [Px,Py]
    % x,y,z scalars

    sigma_z_sqrt_2 = PSF_params.s_ref * sqrt(2 * (1 + (z / PSF_params.z_ref).^2));

    sigma_z_sqrt_2 = reshape(sigma_z_sqrt_2, 1, 1, length(sigma_z_sqrt_2));

    g = 0.25 * pagemtimes( ...
        diff(erf(permute(x_bnd - x', [1 3 2]) ./ sigma_z_sqrt_2)), ...
        diff(erf(permute(y_bnd - y, [3 2 1]) ./ sigma_z_sqrt_2)));
end
