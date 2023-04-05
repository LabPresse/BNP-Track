function U = get_image(fm, Xm, Ym, Zm, ak, x_bnd, y_bnd, dt_exp, F, h, PSF_params)
    % size(fm) = 1 x M
    % size(xm) = K x M
    % size(ym) = K x M
    % size(zm) = K x M
    % size(U) = Px x Py

    if isempty(fm)
        m_idx = 1:size(Xm, 2);
    else
        m_idx = find(fm);
    end

    U = zeros(length(x_bnd) - 1, length(y_bnd) - 1);

    for m = m_idx

        for k = 1:length(ak)
            U = U + ak(k) * get_g_PSF(Xm(k, m), Ym(k, m), Zm(k, m), x_bnd, y_bnd, PSF_params);
        end

    end

    U = dt_exp * (F * diff(x_bnd) * diff(y_bnd) + h * U);
end
