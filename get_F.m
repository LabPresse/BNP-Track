function F = get_F(i,params)

F = params.F_TC + (params.F_IC-params.F_TC) * exp( - i/params.i_sc );
