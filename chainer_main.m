function chain = chainer_main(chain_init,d_length,opts,flag_status,flag_visual)
% to init:
% chain = chainer_main([]   ,  0, opts, true, []  );
% to expand:
% chain = chainer_main(chain,+25, []  , true, true);
% to reduce:
% chain = chainer_main(chain,-10, []  , true, []  );
% to reset record:
% chain = chainer_main(chain, [], []  , true, []  );

tic_id = tic;


% initialize the seed or use seed for last expansion
if isempty(chain_init)
    rng('shuffle');
else
    rng(chain_init.random);
end


% init chain --------------------------------------------------------------
if d_length == 0

    % MCMC
    chain.params = chainer_init_params(opts);
    chain.length = 1;
    chain.stride = 1;
    chain.ledger = nan(0,2); % total wall time without Initialization
    chain.sizeGB = nan;      % current chain memory size
    chain.record = [];       % acceptance rates
    chain.sample = [];

    chain.sample = chainer_init_sample(chain.params,opts);
    
    % history
    chain.i = cast( chain.sample.i    , 'uint64' );
    chain.L = cast( chain.sample.L    , 'double' );
    chain.F = cast( chain.sample.F    , 'single' );
    chain.C = cast( chain.sample.C'   , 'single' ); % reshape
    chain.D = cast( chain.sample.D    , 'single' );
    chain.b = cast( chain.sample.b    , 'logical');
    chain.h = cast( chain.sample.h'   , 'single' ); % reshape
    chain.K = cast( chain.sample.K    , 'uint64' );
    chain.X = cast( chain.sample.X(:)', 'single' ); % reshape
    chain.Y = cast( chain.sample.Y(:)', 'single' ); % reshape
    chain.Z = cast( chain.sample.Z(:)', 'single' ); % reshape

    
    if flag_status
        disp('CHAINER: chain initiated')
    end

    

% expand chain ------------------------------------------------------------
elseif d_length > 0

    chain.params = chain_init.params;
    chain.length = chain_init.length + d_length;
    chain.stride = chain_init.stride;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;
    
    chain.i = [ chain_init.i; zeros( d_length, 1                             , 'like', chain_init.i )];
    chain.L = [ chain_init.L;   nan( d_length, get_log_probs                 , 'like', chain_init.L )];
    chain.F = [ chain_init.F;   nan( d_length, 1                             , 'like', chain_init.F )];
    chain.C = [ chain_init.C;   nan( d_length, chain.params.N                , 'like', chain_init.C )];
    chain.D = [ chain_init.D;   nan( d_length, 1                             , 'like', chain_init.D )];
    chain.b = [ chain_init.b; false( d_length, chain.params.M                , 'like', chain_init.b )];
    chain.h = [ chain_init.h;   nan( d_length, chain.params.N                , 'like', chain_init.h )];
    chain.K = [ chain_init.K; zeros( d_length, chain.params.M                , 'like', chain_init.K )];
    chain.X = [ chain_init.X;   nan( d_length, chain.params.M*chain.params.N , 'like', chain_init.X )];
    chain.Y = [ chain_init.Y;   nan( d_length, chain.params.M*chain.params.N , 'like', chain_init.Y )];
    chain.Z = [ chain_init.Z;   nan( d_length, chain.params.M*chain.params.N , 'like', chain_init.Z )];
    
    
    if flag_visual
        Gim = chainer_visualize([],chain);
    end
    
    %---------------------------- expand chain
    r = chain_init.length+1;
    while r <= chain.length
        
        chain.sample = sampler_update(chain.sample,chain.params);
        
        if mod(chain.sample.i,chain.stride) == 0
            
            chain.i(r  ) =     chain.sample.i;
            chain.L(r,:) =     chain.sample.L;
            chain.F(r  ) =     chain.sample.F;
            chain.C(r,:) =     chain.sample.C;
            chain.D(r  ) =     chain.sample.D;
            chain.b(r,:) =     chain.sample.b;
            chain.h(r,:) =     chain.sample.h;
            chain.K(r,:) =     chain.sample.K;
            chain.X(r,:) =     chain.sample.X(:)'; 
            chain.Y(r,:) =     chain.sample.Y(:)';
            chain.Z(r,:) =     chain.sample.Z(:)';

            
            if flag_visual
                chainer_visualize(Gim,chain);
            end
            
            if flag_status
                disp([  'i = ', num2str(chain.sample.i,'%d'), ...
                     ' - B = ', num2str(sum(double(chain.sample.b)),'%d'), ...
                     ' - acc = ', ...
                                num2str( chain.sample.rec(1,:)./chain.sample.rec(2,:)  * 100 ,'%#6.2f') , ' %', ...
                     ])
            end
            
            r = r+1;
        end
    end    

    if flag_status
        disp('CHAINER: chain expanded')
    end


% reduce chain ------------------------------------------------------------
elseif d_length < 0

    d_length = min(-d_length,chain_init.length);

    chain.params = chain_init.params;
    chain.length = d_length;
    chain.stride = nan;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;
    
    ind = mod(chain_init.length,d_length)+(floor(chain_init.length/d_length)*(1:d_length));

    chain.i = chain_init.i(ind  );
    chain.L = chain_init.L(ind,:);
    chain.F = chain_init.F(ind  );
    chain.C = chain_init.C(ind,:);
    chain.D = chain_init.D(ind  );
    chain.h = chain_init.h(ind,:);
    chain.b = chain_init.b(ind,:);
    chain.K = chain_init.K(ind,:);
    chain.X = chain_init.X(ind,:);
    chain.Y = chain_init.Y(ind,:);
    chain.Z = chain_init.Z(ind,:);

    
    chain.stride = double(chain.i(2)-chain.i(1));
    
    
    if flag_status
        disp('CHAINER: chain reduced')
    end
    
    
% reset chain -------------------------------------------------------------
elseif isempty(d_length)

    chain = chain_init;
    
    chain.record = [chain.record; [chain.sample.i,chain.sample.rec(1,:)./chain.sample.rec(2,:)] ];

    chain.sample.rec(1,:) = 0;
    chain.sample.rec(2,:) = realmin;

    
    if flag_status
        disp('CHAINER: chain reset')
    end
    
end


% store the seed for future expansion
chain.random = rng();



%% book-keeping
chain.sizeGB = get_sizeGB(chain);               % mem size

% ledger
wall_time = toc(tic_id);
chain.ledger = [chain.ledger; double(chain.i(end)), wall_time];

if flag_status
    disp(['( wall time = ',num2str(wall_time),' s, overall wall time = ',num2str(sum(chain.ledger(:,2))),' s )'])
end


end





%% auxiliary functions

function sizeGB = get_sizeGB(chain)
    sizeGB = whos( inputname(1) );
    sizeGB = sizeGB.bytes/1024^3;
end
