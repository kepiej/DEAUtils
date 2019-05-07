function [h,p,Tn,exitflag,bw] = litest2009(X,Y,varargin)
%   Nonparametric similarity of distribution test of two unknown density 
%   functions as described in:
%
%   Qi Li, Esfandiar Maasoumi & Jeffrey S. Racine (2009) "A nonparametric test for equality of distributions with
%   mixed categorical and continuous data", Journal of Econometrics, 148,
%   186-200, DOI: http://dx.doi.org/10.1016/j.jeconom.2008.10.007
%
%   The unknown densities are estimated using kernel density estimation and
%   p-values are obtained through a bootstrapping procedure by resampling from a pooled sub-sample.
%
%   Author: Pieter Jan Kerstens, 2017
%
%   [h,p,Tn,exitflag,bw] = litest2009(X,Y,alpha,nboot)
%       X,Y: (n1 x p) and (n2 x p) matrices containing the two samples of the two unknown
%       densities
%       alpha: significance level (optional, default = 0.05)
%       nboot: number of bootstrap iterations (optional, default = 2000)
%
%       h: h indicates the result of the hypothesis test:
%           h = 0 => Do not reject the null hypothesis at the (100*alpha)% significance level.
%           h = 1 => Reject the null hypothesis at the (100*alpha)% significance level.
%       p: p-value obtained by bootstrapping (see Li et al.(2009) for details)
%       Tn: test statistic
%       exitflag: exitflag from the bandwidth optimization
%       bw: the kernel bandwidth
%
%   See also: smootheffscorebeforelitest


    defopt = {0.05,2000};
    defopt(1:length(varargin)) = varargin;
    [alpha,nboot] = defopt{:};

    assert(size(X,2) == size(Y,2),'X and Y must have the same number of columns p!');
        
    % Kernel function
    kernelft = @(t) exp(-(t.^2)/2)./sqrt(2*pi);
    % Two-fold convolution kernel
    kernelconvft = @(t) exp(-(t.^2)/4)./sqrt(4*pi);
    
    % Function handle to kernel estimator
    kernelest = @(U,V,t) kernelft((repmat(U,1,size(V,1)) - repmat(V,1,size(U,1))')./t)./t;
    kernelconvest = @(U,V,t) kernelconvft((repmat(U,1,size(V,1)) - repmat(V,1,size(U,1))')./t)./t;

    % Pooled sample
    Z = [X;Y];
    N = size(Z,1);
    
    % Least squares leave-one-out cross-validation
    function res = lscv(t)
        temp = prodkernelest(kernelest,Z,Z,t);
        temp(logical(eye(size(temp)))) = 0;
        
        res = (sum(sum(prodkernelest(kernelconvest,Z,Z,t)))/(N^2)) - (2/(N*(N-1)))*sum(sum(temp));
    end

    % Rule-of-thumb bandwidth
    starth = 0.9.*min([std(Z); iqr(Z)/1.349].*N^(-1/5),[],1);
    if(any(starth <= 1e-12))
        starth(starth <= 1e-12) = 1e-3;
    end
    
    % Cross-validation
    opt = optimset('Display','off');
    [bw,~,exitflag] = fminsearch(@lscv,starth,opt);
    if(exitflag ~= 1)
        % Try again with limited search
        [bw,~,exitflag] = fminbnd(@lscv,0.1.*starth,10.*starth);
        if(exitflag < 1)
            warning('litest2009:bw','Cross-validation did not find optimal bandwidth bw! Using rule-of-thumb instead!');
            bw = starth;
        end
    end
    
    % Compute test statistic
    Tn = litesthelper(X,Y,@(U,V) prodkernelest(kernelest,U,V,bw),bw);
    
    % Apply bootstrap procedure on pooled sample Z to approximate
    % empirical distribution of test statistic Tn
    Tnboot = NaN(1,nboot);
    parfor nb=1:nboot
        Tnboot(nb) = litesthelper(datasample(Z,size(X,1)),datasample(Z,size(Y,1)),@(U,V) prodkernelest(kernelest,U,V,bw),bw);
    end
    
    % Determine p-value from the empirical bootstrap distribution
    p = sum(Tnboot > Tn)/nboot;
    h = (p <= alpha);
end

function Tn = litesthelper(X,Y,kernelest,h)
    n1 = size(X,1);
    n2 = size(Y,1);
    
    % Compute kernel estimates and set diagonal entries to 0.
    % This is equivalent to removing the center term in the double sum.
    KXX = kernelest(X,X);
    KXX(logical(eye(size(KXX)))) = 0;
    KYY = kernelest(Y,Y);
    KYY(logical(eye(size(KYY)))) = 0;
    KXY = kernelest(X,Y);
    KYX = kernelest(Y,X);
    
    % Test statistic
    I = (sum(sum(KXX))/(n1*(n1-1))) + (sum(sum(KYY))/(n2*(n2-1))) - ((sum(sum(KXY)) + sum(sum(KYX)))/(n1*n2));
    
    % Estimate of variance
    Omega = 2*n1*n2*prod(h)*((sum(sum(KXX.^2))/((n1^2)*((n1-1)^2))) + (sum(sum(KYY.^2))/((n2^2)*((n2-1)^2))) + ((sum(sum(KXY.^2)) + sum(sum(KYX.^2)))/((n1^2)*(n2^2))));
    
    % Rescaled test statistic which has an asymptotic normal distribution
    Tn = (sqrt(n1*n2*prod(h))*I)/sqrt(Omega);
end

% Compute product kernel
function est = prodkernelest(kernelest,U,V,t)
    est = kernelest(U(:,1),V(:,1),t(1));
    for k=2:length(t)
        est = est.*kernelest(U(:,k),V(:,k),t(k));
    end
end

% Weights for the dependent wild bootstrap (DWB)
function W = DWBweights(n,ln)
    eta = randn(n,1);
    W = zeros(n,n);
    W(1,1) = (exp(-1/ln).*randn) + sqrt(1-exp(-2/ln))*eta(1);
    for k=2:n
        W(k,1:k) = (exp(-1/ln).*W(k-1,1:k)) + sqrt(1-exp(-2/ln))*eta(k);
    end
end