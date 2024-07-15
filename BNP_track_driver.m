clear

%% load data
load('demo_data.mat','w')

%% init chain
opts.units.image  = 'ADU';
opts.units.time   = 's';
opts.units.length = '\mum';

opts.t_ded = 0.01; % dead time

opts.t_exp = 0.10*ones(1,size(w,3)); % exposure periods for each frame

opts.x_bnd = 0.133*(0:size(w,1)); % pixel edges for all pixels
opts.y_bnd = 0.133*(0:size(w,2)); % pixel edges for all pixels

opts.w_cnt = w; % image data

opts.wM    = 200;  % pixel offset
opts.wV    = 16;   % pixel variance
opts.wG    = 0.34; % overall gain
opts.wF    = 2;    % excess noise factor

opts.lambda = 0.565; % imaged wavelength

opts.nRI    = 1.51; % immersion fluid
opts.nNA    = 1.49; % numerical aperture


%% init chain
chain = chainer_main([],0,opts,true,[]);


%% expand chain
chain = chainer_main(chain,500,[],true,false); % expands by adding 500 new samples and reports progress only in command window

