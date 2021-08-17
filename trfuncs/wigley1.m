function [sss,epps]=wigley1(r,NN)
% wigley1:  subsample signal strength and expressed population signal
% [sss,epps]=wigley1(r,NN);
% Last revised 2006-2-24
%
% Subsample signal strength (SSS) and expressed population signal (EPS).  Depending on
% form of r and NN, returns SSS and EPS as sequences paired with given
% constant r and sample sizes from 1 to NN, or as time series (see notes)
%
%*** IN ARGS 
%
% r (1 x 1)r or (? x 1)r -- effective chronolgy signal; a scalar or a time series (Notes)
% NN (1 x 1)i, or (? x 1)i -- maximum number of series at any time, or time series of number of trees (Notes)
%
%*** OUT ARGS 
%
% sss (NN x 1)r or (? x 1)r subsample signal strength as a function of sample size (NN x 1), or as a function of
%       time (? x 1)r
% epps (NN x 1)r or (? x 1)r expressed population signal as a function of sample size (NN x 1), or as a function of
%       time (? x 1)r
%
%
% REFERENCES
%
% Briffa K. and Jones P. D. (1990) Basic chronology statistics and assessment. In Methods of Dendrochronology,
% Applications in the Environmental Sciences (ed E. R. Cook and L. A. Kairiukstis), pp. 136-152. Kluwer Academic Publishers.
%
% Osborn T. J., Briffa K. R. and Jones P. D. (1997) Adjusting variance for sample-size in tree-ring chronologies
% and other regional mean timeseries. Dendrochronologia 15, 89-99.
%
% Wigley T. M. L., Briffa K. R. and Jones P. D. (1984) On the average value of correlated time series, with
% applications in dendroclimatology and hydrometeorology. Journal of Climate and Applied Meteorology 23, 201-213.
%% Wigley et al. (1984)
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES  
%
% In tree-ring application, sample size is number of trees, not
% number of cores.
%
% Variable name epps used instead of eps because eps is a built-in variable
% storing the machine precision
%
% Form of r and NN.  If r and NN are scalars, sss and eps are col vectors of
% length 1:NN, giving, say, sss at 1 tree, 2 trees,..., NN trees.  If r and
% NN are column vectors,r is assumed to be a time series of time-varying
% between-tree correlation, and NN is assumed to be a corresponding time
% series of number of trees in the chronology.
%
% The input effective chronology signal, r, can be computed from the time series and its sample-size
% information by function rbareff.m


%--- EPS : eqn 3.43, Wigley & Jones 1990, p. 146

if all(size(r)==1); % if r a scalar
    if ~all(size(NN)==1);
        error('r is scalar, but NN is not');
    end;
    n = (1:NN)';
    nmax=NN; % maximum sample size)
    top = n*r;
    bottom =  1+ (n-1)*r;
    kwhat='scalar';
else;
    if size(r,2)~=1 | size(NN,2)~=1;
        error('col size of r and NN not 1');
    end;
    if size(r,1) ~= size(NN,1);
        error('r and NN not same row size');
    end;
    n=NN; % number of trees 
    nmax = max(n); % maximum sample size
    top = n .* r;
    bottom =  1+ (n-1) .* r;
    kwhat='vector';
end;
epps = top ./ bottom;
epps(epps<0)=0;



%--- SSS

N = repmat(nmax,length(n),1);

switch kwhat;
    case 'scalar';
        top = n .* (1 + (N-1)*r);
        bottom = N .* (1 + (n-1)*r);
    case 'vector';
        top = n .* (1 + (N-1) .* r);
        bottom = N .* (1 + (n-1) .* r);
    otherwise;
end;
sss = top ./ bottom;

% End of file
