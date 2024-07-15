function log_L = get_log_like(u_cnt,F,params,n)

if isempty(n)
    log_L = - 0.5 * sum( (params.dW_cnt(:)-u_cnt(:) ).^2 ./ ( params.wV/params.wG^2 + F * params.wF*u_cnt(:) ) );
else
    log_L = - 0.5 * sum( (params.dW_cnt(:,:,n)-u_cnt).^2 ./ ( params.wV/params.wG^2 + F * params.wF*u_cnt    ) , [1 2]);
end
