function log_p = log_inv_gam_pdf(x, a, b)
    log_p = a*log(b) - gammaln(a) - (a+1)*log(x)-b./x;
end
