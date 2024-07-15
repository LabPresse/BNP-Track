function D = get_D(X,Y,Z,b,params)


% D|active
m = find(b);
D = params.D_prior_B + 0.25*sum( [diff(X(:,m)),diff(Y(:,m)),diff(Z(:,m))].^2 ./ params.dt,[1 2]);
D = D/randg(params.D_prior_A+1.5*(params.N-1)*length(m));


if nargout<1
    
    Dx = params.D_prior_B + 0.25*sum( diff(X(:,m)).^2 ./ params.dt );
    Dy = params.D_prior_B + 0.25*sum( diff(Y(:,m)).^2 ./ params.dt );
    Dz = params.D_prior_B + 0.25*sum( diff(Z(:,m)).^2 ./ params.dt );

    Dx = Dx / (params.D_prior_A+0.5*(params.N-1)-1);
    Dy = Dy / (params.D_prior_A+0.5*(params.N-1)-1);
    Dz = Dz / (params.D_prior_A+0.5*(params.N-1)-1);

    disp(['Overall: D = ',num2str(D )])
    disp( ' ' )
    disp([' x only: D = ',num2str(Dx)])
    disp([' y only: D = ',num2str(Dy)])
    disp([' z only: D = ',num2str(Dz)])

    clear D
end