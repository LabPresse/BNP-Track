function w = tiff2mat(tifffile)

tiffinfo = imfinfo(tifffile);
nframes = length(tiffinfo);

w = zeros(tiffinfo(1).Width,tiffinfo(1).Height,nframes);

for i = 1:nframes
    w(:,:,i) = fliplr(imread("test.tiff", i)');
end

end