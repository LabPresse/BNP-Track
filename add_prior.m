function add_prior(ax,tag,p1,p2)

% RV
x_temp = linspace(min(get(ax,'Ylim')),max(get(ax,'Ylim')));

% PDF
switch tag
    
    case 'norm' % p1 = mu, p2 = sg
        p_temp = eps + normpdf(x_temp,p1,p2);
    
    case 'symnorm' % p1 = mu, p2 = sg
        p_temp = eps + 0.5*normpdf(x_temp,+p1,p2) + 0.5*normpdf(x_temp,-p1,p2);
    
    case 'foldnorm' % p1 = mu, p2 = sg
        p_temp = eps + normpdf(x_temp,+p1,p2) + normpdf(x_temp,-p1,p2);

    case 'gamma' % p1 = phi, p2 = psi
        p_temp = eps + (p1/p2)^p1/gamma(p1)*x_temp.^(p1-1).*exp(-x_temp*p1/p2);

    case 'beta' % p1 = A, p2 = B
        p_temp = eps + betapdf(x_temp,p1,p2);
    
    case 'invgamma' % p1 = phi, p2 = chi
        p_temp = eps + (p1*p2)^p1/gamma(p1)*x_temp.^(-p1-1).*exp(-p1*p2./x_temp);
    
    case 'symgamma' % p1 = A, p2 = B
        p_temp = eps + 0.5*(p1/p2)^p1/gamma(p1)*abs(x_temp).^(p1-1).*exp(-abs(x_temp)*p1/p2);
    
    case 'besselK' % p1(1:2) = phi(1:2), p2 = psi(1:2)
        p_temp = eps + 2/(gamma(p1(1))*gamma(p1(2)))*(p1(1)*p1(2)/(p2(1)*p2(2))*x_temp).^(0.5*(p1(1)+p1(2))).*besselk(p1(1)-p1(2),2*sqrt(p1(1)*p1(2)/(p2(1)*p2(2))*x_temp))./x_temp;

    otherwise
        warning('Unknown tag (prior visualizer)')
end

% print
p_temp = min(get(ax,'Xlim')) + 0.15*(max(get(ax,'Xlim'))-min(get(ax,'Xlim')))*p_temp/max(p_temp);
line(p_temp,x_temp,'color','k','linestyle','--','parent',ax);


