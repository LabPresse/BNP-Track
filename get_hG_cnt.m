function hG_cnt = get_hG_cnt(hP_h_t_exp,x,y,z,x_bnd,y_bnd,s_ref,z_ref)
% hG_cnt = integral of emitter PSF over pixels
% x,y,z are scalars

sigma_sqrt_2 = sqrt(2)*s_ref*sqrt( 1 + (z/z_ref)^2 );
sqrt_hpt_2 = 0.5*sqrt( hP_h_t_exp );

hG_cnt = ( sqrt_hpt_2 * diff(erf ((x_bnd-x)/sigma_sqrt_2)) ) ...
       * ( sqrt_hpt_2 * diff(erf ((y_bnd-y)/sigma_sqrt_2)) );
