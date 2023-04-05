function PSF_params = get_PSF_params(NA, n, lambda)

    if nargin < 1
        NA = 1.45; % numerical aperture
        n = 1.515; % immersion fluid index
        lambda = 561/1000; % [length] emission wavelength in vacuum
    end

    cosalphasqrt = sqrt(cos(asin(NA / n)));
    cos32alpha = cosalphasqrt^3;
    cos72alpha = cosalphasqrt^7;

    PSF_params.tag = 1; % 'gauss';
    PSF_params.z_ref = lambda / pi / n * ((7 * (1 - cos32alpha)) / (4 - 7 * cos32alpha + 3 * cos72alpha)); % [length] std of PSF along z  (optical axis)
    PSF_params.s_ref = 0.5 * lambda / pi / n * sqrt((7 * (1 - cos32alpha)) / (4 - 7 * cos32alpha + 3 * cos72alpha)); % [length] std of PSF along xy (image plane)

    % PSF_params.tag   = 'gauss';
    % PSF_params.z_ref = 218.76/1000;   % [length] std of PSF along z  (optical axis)
    % PSF_params.s_ref =  77.35/1000;   % [length] std of PSF along xy (image plane)

    if nargout < 1
        disp(PSF_params)
        clear PSF_params
    end
