function g = get_g_PSF(x, y, z, x_bnd, y_bnd, PSF_params)
    % g = int_Ap G(x,y,z;x,y) dx dt
    % size(g) = [Px,Py]
    % x,y,z scalars

    sigma_z_sqrt_2 = PSF_params.s_ref * sqrt(2 * (1 + (z / PSF_params.z_ref) ^ 2));

    g = 0.25 ...
        * diff(erf((x_bnd - x) / sigma_z_sqrt_2)) ...
        * diff(erf((y_bnd - y) / sigma_z_sqrt_2));

end
