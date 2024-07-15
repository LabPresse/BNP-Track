function Gim = chainer_visualize(Gim,chain)

chain_B = sum(chain.b,2);

chain_b = repmat(1:chain.params.M,chain.length,1);
chain_b(~chain.b) = nan;

X = chain.sample.X;
Y = chain.sample.Y;
Z = chain.sample.Z;
idx = find(~chain.sample.b);
X(:,idx) = nan;
Y(:,idx) = nan;
Z(:,idx) = nan;

k_idx = sub2ind([chain.params.N chain.params.M],chain.sample.K,1:chain.params.M);
K = chain.sample.K;
x = X(k_idx);
y = Y(k_idx);
z = Z(k_idx);


%% init
if isempty(Gim)
    
    col_1 = [  0 204   0]/255; % bright
    col_2 = [ 51 102 204]/255; % dark
    col_3 = [255 204   0]/255; % dark
    
    num = 7;
    mum = 5;
    
    chain_i = double(chain.i(1)) + chain.stride*(0:chain.length-1)';
    i_lim = [max(chain_i(1),0.1*chain_i(end)) chain_i(end)+1];


    
    figure(10)
    clf

    
    
    % --- Sample ----------------------------------------------------------
    ax_X = subplot(num,mum,[1 mum*(num-3)+1]);
    ax_Y = subplot(num,mum,[2 mum*(num-3)+2]);
    ax_Z = subplot(num,mum,[3 mum*(num-3)+3]);

    Gim.X = plot(ax_X,X,chain.params.t_mid,'.-','color',col_1);
    Gim.Y = plot(ax_Y,Y,chain.params.t_mid,'.-','color',col_1);
    Gim.Z = plot(ax_Z,Z,chain.params.t_mid,'.-','color',col_1);

    Gim.x = line(ax_X,x,chain.params.t_mid(K),'marker','*','linestyle','none','color',col_3);
    Gim.y = line(ax_Y,y,chain.params.t_mid(K),'marker','*','linestyle','none','color',col_3);
    Gim.z = line(ax_Z,z,chain.params.t_mid(K),'marker','*','linestyle','none','color',col_3);

    ax_X.Box = 'off';
    ax_Y.Box = 'off';
    ax_Z.Box = 'off';
    
    xline(ax_X,chain.params.x_bnd([1 end]),':')
    xline(ax_Y,chain.params.y_bnd([1 end]),':')

    yline(ax_X,chain.params.t_bnd,':')
    yline(ax_Y,chain.params.t_bnd,':')
    yline(ax_Z,chain.params.t_bnd,':')

    xline(ax_X,chain.params.X_prior_min,'--')
    xline(ax_Y,chain.params.Y_prior_min,'--')
    xline(ax_Z,chain.params.Z_prior_min,'--')
    
    xline(ax_X,chain.params.X_prior_max,'--')
    xline(ax_Y,chain.params.Y_prior_max,'--')
    xline(ax_Z,chain.params.Z_prior_max,'--')

    xline(ax_Z,0,':')

    ax_X.XLim = [min(ax_X.XLim(1),min(X(:))) max(ax_X.XLim(2),max(max(X(:))))];
    ax_Y.XLim = [min(ax_Y.XLim(1),min(Y(:))) max(ax_Y.XLim(2),max(max(Y(:))))];
    ax_Z.XLim = [min(ax_Z.XLim(1),min(Z(:))) max(ax_Z.XLim(2),max(max(Z(:))))];

    ax_X.YLim = chain.params.t_bnd([1 end]);
    ax_Y.YLim = chain.params.t_bnd([1 end]);
    ax_Z.YLim = chain.params.t_bnd([1 end]);
    
    clim(ax_X,'manual')
    clim(ax_Y,'manual')
    clim(ax_Z,'manual')
    
    xlabel(ax_X,['X^{1:M} (',chain.params.units.length,')'])
    xlabel(ax_Y,['Y^{1:M} (',chain.params.units.length,')'])
    xlabel(ax_Z,['Z^{1:M} (',chain.params.units.length,')'])

    ylabel(ax_X,['T (',chain.params.units.time,')'])
    
    ax_Y.YTickLabel = [];
    ax_Z.YTickLabel = [];

    title(ax_Y,'MCMC sample')

    

    ax_T = subplot(num,mum,[mum*(num-2)+1 mum*num-2]);

    Gim.T = plot3(ax_T,X,Y,Z,'.-','color',col_1);
    ax_T.Box = 'off';
 
    grid(ax_T,'on')

    xlabel(ax_T,['X^{1:M} (',chain.params.units.length,')'])
    ylabel(ax_T,['Y^{1:M} (',chain.params.units.length,')'])
    zlabel(ax_T,['Z^{1:M} (',chain.params.units.length,')'])

    ax_T.XLim = [min(ax_X.XLim(1),min(X(:))) max(ax_X.XLim(2),max(max(X(:))))];
    ax_T.YLim = [min(ax_Y.XLim(1),min(Y(:))) max(ax_Y.XLim(2),max(max(Y(:))))];
    ax_T.ZLim = [min(ax_Z.XLim(1),min(Z(:))) max(ax_Z.XLim(2),max(max(Z(:))))];

    % draw prior box
    R = [chain.params.X_prior_min, chain.params.X_prior_max, chain.params.X_prior_max, chain.params.X_prior_min,        chain.params.X_prior_min, chain.params.X_prior_max, chain.params.X_prior_max, chain.params.X_prior_min;
         chain.params.Y_prior_min, chain.params.Y_prior_min, chain.params.Y_prior_max, chain.params.Y_prior_max,        chain.params.Y_prior_min, chain.params.Y_prior_min, chain.params.Y_prior_max, chain.params.Y_prior_max;
         chain.params.Z_prior_min, chain.params.Z_prior_min, chain.params.Z_prior_min, chain.params.Z_prior_min,        chain.params.Z_prior_max, chain.params.Z_prior_max, chain.params.Z_prior_max, chain.params.Z_prior_max];
    line(ax_T,[R(1,1) R(1,2) R(1,3) R(1,4) R(1,1)],...
              [R(2,1) R(2,2) R(2,3) R(2,4) R(2,1)],...
              [R(3,1) R(3,2) R(3,3) R(3,4) R(3,1)],...
         'color','k','linestyle',':')
    line(ax_T,[R(1,5) R(1,6) R(1,7) R(1,8) R(1,5)],...
              [R(2,5) R(2,6) R(2,7) R(2,8) R(2,5)],...
              [R(3,5) R(3,6) R(3,7) R(3,8) R(3,5)],...
         'color','k','linestyle',':')
    line(ax_T,[R(1,1) R(1,5)],...
              [R(2,1) R(2,5)],...
              [R(3,1) R(3,5)],...
         'color','k','linestyle',':')
    line(ax_T,[R(1,2) R(1,6)],...
              [R(2,2) R(2,6)],...
              [R(3,2) R(3,6)],...
         'color','k','linestyle',':')
    line(ax_T,[R(1,3) R(1,7)],...
              [R(2,3) R(2,7)],...
              [R(3,3) R(3,7)],...
         'color','k','linestyle',':')
    line(ax_T,[R(1,4) R(1,8)],...
              [R(2,4) R(2,8)],...
              [R(3,4) R(3,8)],...
         'color','k','linestyle',':')
    
    % draw image box
    line(ax_T,[chain.params.x_bnd(1) chain.params.x_bnd(end) chain.params.x_bnd(end) chain.params.x_bnd(  1) chain.params.x_bnd(1)],...
              [chain.params.y_bnd(1) chain.params.y_bnd(  1) chain.params.y_bnd(end) chain.params.y_bnd(end) chain.params.y_bnd(1)],'color','k')





    

    
    % --- MCMC ------------------------------------------------------------
    ax_F = subplot(num,mum,0*mum+[4 1*mum],'YAxisLocation','Right','XLim',i_lim);
    ax_b = subplot(num,mum,1*mum+[4 2*mum],'YAxisLocation','Right','XLim',i_lim);
    ax_C = subplot(num,mum,3*mum+[4 1*mum],'YAxisLocation','Right','XLim',i_lim);
    ax_h = subplot(num,mum,4*mum+[4 1*mum],'YAxisLocation','Right','XLim',i_lim);
    ax_D = subplot(num,mum,5*mum+[4 2*mum],'YAxisLocation','Right','XLim',i_lim);
    
    title(ax_F,['MCMC chain (stride=',num2str(chain.stride),')'])
	xlabel(ax_D, 'MCMC iteration (i)')

    ax_F.YScale = 'log';
    ax_C.YScale = 'log';
    ax_h.YScale = 'log';
    ax_D.YScale = 'log';

    ax_F.YGrid = 'on';
    ax_C.YGrid = 'on';
    ax_h.YGrid = 'on';
    ax_D.YGrid = 'on';

    ax_F.XTickLabel = [];
    ax_b.XTickLabel = [];
    ax_C.XTickLabel = [];
    ax_h.XTickLabel = [];

    ylabel(ax_F, 'F (1)' )
    ylabel(ax_b,{'B, b^{1:M}',[]})
    ylabel(ax_C,['\langleC\rangle (phts/',chain.params.units.time,'/',chain.params.units.length,'^2)'])
    ylabel(ax_h,['\langleh\rangle (phts/',chain.params.units.time,')'])
    ylabel(ax_D,['D (',chain.params.units.length,'^2/',chain.params.units.time,')'])

    ax_b.YLim = [0 chain.params.M+1];

    if ~isempty( chain.ledger )
        xline(ax_F,chain.ledger(end,1));
        xline(ax_b,chain.ledger(end,1));
        xline(ax_C,chain.ledger(end,1));
        xline(ax_h,chain.ledger(end,1));
        xline(ax_D,chain.ledger(end,1));
    end

    Gim.F = line(ax_F,chain_i,      chain.F(:,1),'marker','.','color',col_2);
    Gim.b = line(ax_b,chain_i,     chain_b      ,'marker','.','color',col_1);
    Gim.B = line(ax_b,chain_i,     chain_B      ,'marker','o','color',col_3);
    Gim.C = line(ax_C,chain_i,mean(chain.C,2)   ,'marker','.','color',col_1);
    Gim.h = line(ax_h,chain_i,mean(chain.h,2)   ,'marker','.','color',col_1);
    Gim.D = line(ax_D,chain_i,     chain.D      ,'marker','.','color',col_3);


    if isfield( chain.params,'ground' )
        yline(ax_D,chain.params.ground.D,'color',col_3)
    end
    
end % init





Gim.x.XData = x; Gim.x.YData = chain.params.t_mid(K);
Gim.y.XData = y; Gim.y.YData = chain.params.t_mid(K);
Gim.z.XData = z; Gim.z.YData = chain.params.t_mid(K);

for m = 1:chain.params.M
    % --- sample
    Gim.T(m).XData = X(:,m);
    Gim.T(m).YData = Y(:,m);
    Gim.T(m).ZData = Z(:,m);

    Gim.X(m).XData = X(:,m);
    Gim.Y(m).XData = Y(:,m);
    Gim.Z(m).XData = Z(:,m);
    % --- MCMC
    Gim.b(m).YData = chain_b(:,m);
end
    Gim.F.YData = chain.F;
    Gim.B.YData = chain_B;
    Gim.C.YData = mean(chain.C,2);
    Gim.h.YData = mean(chain.h,2);
    Gim.D.YData = chain.D;



drawnow


end % visualize
