function [A,x,y] = GRIDobj2mat(DEM)

% convert GRIDobj to matrix and coordinate vectors
%
% Syntax
%
%     [dem,X,Y] = GRIDobj2mat(DEM)
%
% Description
%
%     convert GRIDobj to matrix and coordinate vectors 
%
% Input
%
%     DEM       instance of GRIDobj class
% 
% Output
%
%     dem       matrix 
%     x,y       coordinate vectors
%                  
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 27. December, 2012

A = DEM.Z;
[x,y] = refmat2XY(DEM.refmat,DEM.size);
