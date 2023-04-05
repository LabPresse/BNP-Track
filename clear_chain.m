function chain = clear_chain(infile, outfile)
    load(infile);
    chain.length = 1;
    chain.ledger = [];
    chain.record = [];
    chain.i = 1;
    chain.T = 400;

    chain.bm = chain.bm(1, :);
    chain.Xm = chain.Xm(1, :);
    chain.Ym = chain.Ym(1, :);
    chain.Zm = chain.Zm(1, :);
    chain.D = chain.D(1);
    chain.F = chain.F(1);
    chain.h = chain.h(1);
    chain.G = chain.G(1);

    chain.sample.i = 0;
    chain.sample.T = 400;
    chain.sample.bm = chain.bm;
    L = chain.params.N * chain.params.K;
    chain.sample.Xm = double(reshape(chain.Xm, L, chain.params.M));
    chain.sample.Ym = double(reshape(chain.Ym, L, chain.params.M));
    chain.sample.Zm = double(reshape(chain.Zm, L, chain.params.M));
    chain.sample.D = chain.D;
    chain.sample.F = chain.F;
    chain.sample.h = chain.h;
    chain.sample.G = chain.G;
    chain.sample.rec_Fh(1, :) = 0;
    chain.sample.rec_Fh(2, :) = eps;
    chain.sample.rec_XYZ(1, :) = 0;
    chain.sample.rec_XYZ(2, :) = eps;

    chain.params.MH_sc = [[100 100 2] [0.05 0.05 0.05]];
    chain.params.i_skip = 1;

    save_file = outfile;
    clear L outfile
    save(save_file)
end