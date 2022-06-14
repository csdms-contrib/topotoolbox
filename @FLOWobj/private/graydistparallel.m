function I = graydistparallel(I,ix)

SEEDS   = false(size(I));
SEEDS(ix) = true;
regions = regionprops(~isinf(I),{'SubarrayIdx','Image'});

for r = 1:numel(regions)
    temp = SEEDS(regions(r).SubarrayIdx{1},regions(r).SubarrayIdx{2});
    regions(r).SEED = find(temp.*regions(r).Image);
    regions(r).Image = I(regions(r).SubarrayIdx{1},regions(r).SubarrayIdx{2});
end

parfor r = 1:numel(regions)
    regions(r).Image = graydist(regions(r).Image,regions(r).SEED,'q');
end

% Map values back to D
for r = 1:numel(regions)
    temp = I(regions(r).SubarrayIdx{1},regions(r).SubarrayIdx{2});
    temp(~isinf(regions(r).Image)) = regions(r).Image(~isinf(regions(r).Image));
    I(regions(r).SubarrayIdx{1},regions(r).SubarrayIdx{2}) = temp;
end
