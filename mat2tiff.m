function mat2tiff(tifffile, w)

w2 = uint16(w);

imwrite(w2(:, end:-1:1, 1)', tifffile);

for j = 2:size(w2, 3)
    imwrite(w2(:, end:-1:1, j)', tifffile, 'WriteMode', 'append');
end

end
