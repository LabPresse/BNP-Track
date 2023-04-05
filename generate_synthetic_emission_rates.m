function [h, F] = generate_synthetic_emission_rates(dt_exp, dx, dy, PSF_params)

    g_temp = get_g_PSF(0, 0, 0, [-0.5 * dx +0.5 * dx], [-0.5 * dy +0.5 * dy], PSF_params);

    h = 5 * 2 * 0.1 * 100 / dt_exp / g_temp; % emitter photons/[time]
    F = 2 * 0.2 * h * g_temp / (dx * dy); % background photons 1/[area*time]

    %h = h /20;
    %F = F / 3.5;

    disp(['Emitter photons over one pixel over one exposure = ', num2str(h * g_temp * dt_exp)])
    disp(['Backgrd photons over one pixel over one exposure = ', num2str(F * dx * dy * dt_exp)])
