function WI = twi(SCA,G,type)

% calculate different topographic wetness indices
%
% Syntax
%     WI = twi(A,G)
%
% Description
%
%     twi calculates the topographic wetness index of Beven and Kirkby
%     (1979). 
%
% Equation
%
%     WI = log(A/tan(beta))
%
% Input 
%
%     A      specific catchment area
%     G      gradient (as returned by gradient8 as tangent)
%     type   'twi': topographic wetness index (Beven and Kirkby, 1979)
%            'saga': SAGA topographic wetness index
%            'rusle': length-slope factor used in the Revised Universal
%            Soil Loss equation (RUSLE), Moore and Burch (1986). Only for
%            slope lengths < 100 m and slopes < 14° (Wilson and Galant, 
%            2000) 
%
% Output
%
%     WI     wetness index
%
% References:
%
% Stefano, C. D., Ferro, V. & Porto, P., 2000: Length Slope Factors for 
% applying the Revised Universal Soil Loss Equation at Basin Scale in 
% Southern Italy. Journal of Agricultural Engineering Research, 75, 
% 349 - 364. DOI: 10.1006/jaer.1999.0514
%
% Wilson, J. P., Gallant, J. C., 2000: Digital terrain analysis, p. 1-27.
% In J.P. Wilson and J. C. Gallant (ed.). Terrain analysis: Principles and
% applications. John Wiley & Sons, New York.
%
% See also: flowacc, gradient8
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 5. February, 2010


error(nargchk(2, 3, nargin))

if nargin == 2;
    type = 'twi';
else
    twis = {'twi','saga','rusle'};
    if ~ any(strcmpi(type,twis))
        error('unknown type')
    end
end
        
    


% apply threshold to gradient
G = max(G,0.001);

switch lower(type)
    case 'twi'
        WI = log(SCA./G);
        
    case 'saga'
        
        beta = atan(G);
        
        iter = true;
        
        counter = 1;
        while iter && counter <= inf;
            SCAmax   = imdilate(SCA,ones(3)).*((1/15).^(beta.*exp(15.^beta)));
            I        = SCA < SCAmax;
            if any(I(:))
                SCA(I)   = SCAmax(I);
                iter     = true;
                counter  = counter+1;
            else
                iter     = false;
                
            end
        end


        WI = log(SCA./G);
        
    case 'rusle'
        
        % length-slope factor (Moore and Burch)
        n  = 1.3;
        m  = 0.4;        
        WI = (m+1)*(SCA/22.13).^m .* (sin(atan(G))/0.0896).^n;
        slope_thres = tand(14);
        WI(SCA<10 & G<slope_thres) = nan;

end
        