function qc = quadratcount(P,varargin)

%QUADRATCOUNT Quadrat count and chi2 test
%
% Syntax
%
%     qc = quadratcount(P)
%     qc = quadratcount(P,pn,pv,...)
%
% Description
%
%     quadratcount partitions the stream network and calculates the number
%     of points in each partition. 
% 
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
%     Parameter name/value pairs
%
%     'seglength'   distance of partitions. 
%     'chi2test'    {true} or false. Use chi squared test to test whether
%                   the point pattern is randomly distributed. 
%     'voronoi'     {true} or false. See PPS/histogram for details.
%
% Output arguments
%
%     qc      structure array with following fields
%             .labelnal  node-attribute list with partitions (quadrants)
%             .qlength   length of each quadrant
%             .qcount    number of points in each quadrant
%             .qint      density of points in each quadrant
%             .qcountexp expected number of counts in each quadrant under
%                        the assmumption of complete spatial randomness
%                        (CSR).
%             .X2hat     estimated chi2 value
%             .X2theo    theoretical chi2 value
%             .X2p       p-value (probability of erroneously rejecting the
%                        null-hypothesis of CSR.
%
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

% Parse inputs
p = inputParser;
p.FunctionName = 'PPS/quadratcount';
addParameter(p,'seglength',tlength(P)/npoints(P)*10);
addParameter(p,'chi2test',true);
addParameter(p,'voronoi',true);
parse(p,varargin{:});

% distance between nodes
d = zeros(size(P.S.x));
d(P.S.ix) = nodedistance(P);

% calculate histogram
[~,lab,cc] = histogram(P,'voronoi',p.Results.voronoi,...
                          'seglength',p.Results.seglength,...
                          'normalization','counts');
nrlab = max(lab);

qc.labelnal = lab;
qc.qlength  = accumarray(lab,d,[nrlab 1]);
qc.qcount   = cc;
qc.qint     = qc.qcount./qc.qlength;

if p.Results.chi2test
    totlength = tlength(P);
    
    qc.qcountexp = zeros(size(qc.qcount));
    for r = 1:numel(qc.qlength)
        PD = makedist('binomial','N',numel(P.PP),'p',qc.qlength(r)./totlength);
        qc.qcountexp(r) = mean(PD);
    end
    
    % There might be segments with zero length
    I = qc.qlength ~= 0;
    n = nnz(I);
    qc.X2hat  = sum(((qc.qcountexp(I)-qc.qcount(I)).^2)./qc.qcountexp(I));
    qc.X2theo = chi2inv(0.99,n-1);
    qc.X2p    = 1-chi2cdf(qc.X2hat,n-1); 
end

    