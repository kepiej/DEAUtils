function TES = smootheffscorebeforelitest(TE,ninputs,noutputs,varargin)
%   Algorithm II described in:
%
%    Léopold Simar & Valentin Zelenyuk (2006) On Testing Equality of
%    Distributions of Technical Efficiency Scores, Econometric Reviews, 25:4, 497-522, DOI:10.1080/07474930600972582
%
%   which "smooths" the efficiency scores of efficient DMUs away from the boundary.
%   This should be used as a preprocessing step when using DEA efficiency scores before using the Li test.
%
%   Author: Pieter Jan Kerstens, 2016
%
%   TES = smootheffscorebeforelitest(TE,ninputs,noutputs,alpha,effval)
%       TE: (k x 1) vector of DEA efficiency scores
%       ninputs: number of inputs
%       noutputs: number of outputs
%       alpha: significance level of the Li test (optional, default = 5%)
%       effval: value that indicates an efficient DMU (optional, default = 1.0)
%
%       TES: (k x 1) vector of "smoothed" DEA efficiency scores
%
%   See also: litest


    defopt = {5,1.0};
    defopt(1:length(varargin)) = varargin;
    [alpha,effval] = defopt{:};

    k = size(TE,1);
    
    % Find DMUs with efficiency score = effval
    ind = find(abs(effval - TE) <= 1e-14);
    
    % Smooth efficient DMUs
    TES = TE;
    TES(ind) = TE(ind) + min(k^(-2/(ninputs+noutputs+1)),alpha-1).*rand(length(ind),1);
end