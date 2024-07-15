function sample = chainer_init_sample(params,opts)


%% counter
sample.i = 0;
sample.F = get_F(sample.i,params);


%% loads
sample.b = rand(1,params.M) < 1/(1+exp(-params.b_prior_log_p1_m_log_p0));


%% emission rates
sample.h = params.h_prior_REF./params.h_prior_A.*randg(params.h_prior_A,params.N,1);


%% atoms and dynamics
sample.K = nan(       1,params.M);
sample.X = nan(params.N,params.M);
sample.Y = nan(params.N,params.M);
sample.Z = nan(params.N,params.M);

[sample.K,sample.X,sample.Y,sample.Z] = sampler_update_DDD(sample.K,sample.X,sample.Y,sample.Z,false(1,params.M),params);


%% background
sample.C = params.C_prior_REF./params.C_prior_A.*randg(params.C_prior_A,params.N,1);


%% book-keeping
sample.L = get_log_probs(sample.F,sample.C,sample.b,sample.h,sample.K,sample.X,sample.Y,sample.Z,params);
sample.D = get_D(sample.X,sample.Y,sample.Z,sample.b,params);

sample.rec = repmat([0;realmin;0],1,9);


