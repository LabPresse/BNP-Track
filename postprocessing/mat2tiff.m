function mat2tiff(datafilename, tiffname)

    load(datafilename, "chain");
    w_cnt = chain.params.w_cnt;
    newmap = gray(256);
    imwrite(w_cnt(:, end:-1:1, 1)', newmap, tiffname);

    for j = 2:size(w_cnt, 3)
        imwrite(w_cnt(:, end:-1:1, j)', newmap, tiffname, 'WriteMode', 'append');
    end

end
