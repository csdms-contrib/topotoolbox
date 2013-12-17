function writeinmatrix(app,dem,ix)

[row,col] = ind2sub(app.DEM.size,ix); 
app.DEM.Z(min(row):max(row),min(col):max(col)) = dem.Z;
end