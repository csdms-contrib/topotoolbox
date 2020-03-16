function varargout=neighbourIndices(A,varargin)
% function part of HyLands: Hybrid Landscape evolution model 
%
% See also: HYLANDS_set, HYLANDS
% 
% =========================================================================
% Papers to cite when using HyLands:
%
% * HyLands: Campforts B., Shobe M.C., et al. : HyLands 1.0: a Hybrid
% Landscape evolution model to simulate the impact of landslides and
% landslide-derived sediment on landscape evolution. Discussion paper in
% Geoscientific Model Development,
% https://geoscientific-model-development.net
%
% * TTLEM: Campforts, B., Schwanghart, W., & Govers, G. (2017). 
% Accurate simulation of transient landscape evolution
% by eliminating numerical diffusion: the TTLEM 1.0 model. Earth Surface
% Dynamics, 5(1), 47–66. https://doi.org/10.5194/esurf-5-47-2017
%
% =========================================================================
%
% Author:   Benjamin Campforts (benjamin.campforts@gfz-potsdam.de)
%
% Date:     15. March, 2020
if ~isempty(varargin)
    validStrings = {'M','ind','M&ind'};
    try
        outputStyle= validatestring(varargin{1},validStrings);
    catch
        error('Valid output types for neighbourIndices are:\n%s',' ''M'', ''ind'' or ''M&ind''');
    end
end
% Structure= [1;2;3;4;5;6;7;8]
% --------- %
% -8--4--5- %
% --------- %
% -3--X--1- %
% --------- %
% -7--2--6- %
% --------- %
% Example
% A=[1 2 3 ; 4 5 6; 7 8 9]
% nb=(neighbourIndices(A));
% neighbours=nan(size(nb));
% neighbours(nb~=0)=A(nb(nb~=0))

nb_row=size(A,1);
nb_col=size(A,2);
nb_tot=numel(A);

neighbM=zeros(numel(A),8);
neighbM(:,1)=(1:nb_tot)'+nb_row;
neighbM(:,2)=(1:nb_tot)'+1;
neighbM(:,3)=(1:nb_tot)'-nb_row;
neighbM(:,4)=(1:nb_tot)'-1;
neighbM(:,5)=(1:nb_tot)'+nb_row-1;
neighbM(:,6)=(1:nb_tot)'+nb_row+1;
neighbM(:,7)=(1:nb_tot)'-nb_row+1;
neighbM(:,8)=(1:nb_tot)'-nb_row-1;



neighbM(neighbM<=0)=0;
neighbM(neighbM>nb_tot)=0;
% Remove periodical pointers in top and bottom rows:
%Top
neighbM(1:nb_row:nb_row*(nb_col-1),[4 5 8])=0;
%Not
neighbM(nb_row:nb_row:nb_row*nb_col,[2 6 7])=0;
neighbM=int64(neighbM);
if isempty(varargin)
    varargout{1}=neighbM;
else
    switch outputStyle
        case 'M'
            varargout{1}=neighbM;
        case 'ind'
            [ic, icd] = calcIcIcd(neighbM);
            varargout{1}=ic;
            varargout{2}=icd;
        case 'M&ind'
            [ic, icd] = calcIcIcd(neighbM);
            varargout{1}=neighbM;
            varargout{2}=ic;
            varargout{3}=icd;
    end
end
    function [IC, ICD] = calcIcIcd(neighbM)
        IC=int64(repmat(1:size(neighbM,1),1,8));
        IC=IC(neighbM~=0)';
        ICD=neighbM(neighbM~=0);
    end
end