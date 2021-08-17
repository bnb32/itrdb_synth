function Result=sitechron1(X,yrX,minlap,kopt,nms)
% sitechron: Compute site chronology from core indices
% function Result=sitechron1(X,yrX,minlap,kopt,nms);
% Last Revised 2006-7-25
%
% Compute site chronology from core indices.
% Optionally stabilizes variance to compensate for time-varying sample size.
% First needed by rwlprep1.m as part of sequence in flowrec reconstructions. 
%
%*** INPUT
%
% X (mX x nX)r core indices (Notes)n-internal NaNs (Notes)
% yrx (mX x 1)i year vector for X
% minlap (1 x 1)i  minimum allowable sample size for an acceptable correlation for
%   mean between-series r computation (may be [] if kopt(1)==2).  
% kopt (1 x 2)i options
%   kopt(1) variance stabilization for time varying sample size
%       ==1 Yes (Notes)
%       ==2 No
%   kopt(2) variance equalization of core indices before making site chronology
%       ==1 Yes (Notes)
%       ==2 No
% nms {? x 1}s or []:  core ids, or [] if general time series use (Notes)
% 
%*** OUTPUT
%
% Result -- structure of results
%   .y (my x 1)r  site chronology, after optional variance stabilization
%   .yry (my x 1)r year vector for y
%   .y1 (my x 1)r  site chronology, before optional variance stabilization
%   .yry1 (my x 1)r year vector for y1
%   .ylo (my x 1)r lowest core index in each year
%   .yhi (my x 1)r highest core index in each year
%   .ys (my x 1)r standard deviation of core index in each year 
%   .SSSt (my x 1)r time dependent subsample signal strength
%   .EPSt (my x 1)r time dependent expressed population signal
%   .rbar (my x 1)r time dependent mean between series correlation
%   .ncores (my x 1)i number of series (cores) averaged each year
%   .ntrees (my x 1)i number of trees each year ([] if nms is [])
%   .rbs (1 x 1)r  mean between-series correlation computed from all overlapping pairs (overlap at least
%       minlap)
%   .SSS (n1 x 1) subsample signal strength from scalar correlation .rbs, as a function of sample size 1 to n1,
%       where n1 is the maximum number of cores in any year
%   .EPS (n1 x 1) expressed population signal, computed from .rbs, like SSS
%   .f (my x 1)r  scaling factor (Osborne and Briffa, 1992) used for variance stabilization
%           all ones  if kopt(2)==2
% 
%function [r,yrr,n1,n2,n3,eflag]=corrtsm1(X,yrX,minover,kopt,mwindow,tshift)

% [r,yrr,n1,n2,n3,eflag]=corrtsm1(X,yrX,minover,kopt,mwindow,tshift);
% Last revised 12-15-03

%
%*** REFERENCES --- none
%
%*** TOOLBOXES NEEDED -- none
%
%*** UW FUNCTIONS CALLED 
%
% corrtsm1: mean between-series correlation for a time series matrix
% wigley1:  subsample signal strength and expressed population signal

%
%*** NOTES
%
% X:  May be standard or residual core indices.  Missing values allowed to be interspersed. Site index
%   is average over number of cores available in any year. Special case of X with just one column (a single
%   core) returns time series with same coverage as X after trimming any leading and trailing NaNs, and 
%   scalars for other output. Imbedded NaNs in X result in NaNs or zero for the corresponding years in the
%   time series output, depending on context. Specifically, the output for single-column X is 
%       Result.y=Result.y1=Result.ylo=Result.yhi=X;  the core index is returned as the site chronology, etc
%       Result.yry=yrX, Result.yry1=yrX
%       Result.ys=time series of zeros for the standard deviation in each year
%       Result.ntrees, .ncores, .SSSt, .f = time series of ones
%       Result.EPSt, .rbar = time series of NaNs
%       Result.SSS = 1 (scalar)
%       Result.EPS = NaN
%
%
% Variance stabilization.  If kopt(1)==1, variance is stabilized by method described by Osborne and Briffa
% (1992).  The time-varying correlation method is used.  Correlations are first computed between all pairs of
% series for all overlaps longer than minlat years. Mean between series correlation is average of correlations
% for all pairs of series available in a given year.  See corrtsm1.m for between series correlation method. 
%
% Variance equalization.  if kopt(2)==1, core indices are scaled to have same standard deviation before
% averaging to form site chronology.  The standard deviation of each core index is set to the median standard
% deviation of all core indices.  This is an important step in variance stabilization by the methods
% recommended by Osborne and Briffa (1992), as it is assumed all core indices have equal variance. 
%
% Rescaling of final chronology index.  If variance stabilization is on (kopt(1)==1), the variance-stablilized
% chronology is re-scaled to exactly a mean of 1.0 and a variance equal to the variance of the chronology
% before variance-stabilization.  Rescaling helps in removing differences of level and spread in comparing
% versions of the chronology variance rescaled and not. 
%
% nms:  if [], the statistics SSS, EPS, rbar, rbs are based on correlations between series.  If nms not [],
% the correlations used are only those between trees.   ntrees is returned [] if nms is [].


%--- CHECK INPUT

% X
[mX,nX]=size(X);
if (size(yrX,1) ~= mX) | (size(yrX,2) ~=1);
    error('yrX must be col vector same row size as X');
end
if isempty(X);
    error(['X is empty -- no cores in set']);
end


% kopt
if size(kopt,1) ~= 1;
    error('kopt must be scalar or row vector');
end;
if size(kopt,2) ~=2; 
    error('kopt must be length 2 ');
end
if ~any(kopt(1) == [1 2]);
    error('kopt(1) must be 1 or 2');
end;
if kopt(1)==1;
    if isempty(minlap) | (minlap > mX);
        error(['minlap cannot be empty or larger than row size of X']);
    else;
        minover=minlap; % call to rbareff used the name minover
    end;
end
if ~any(kopt(2) == [1 2]);
    error('kopt(2) must be 1 or 2');
end;


%--- TRIM OFF ANY ALL-NAN TRAILING OR LEADING ROWS

if nX~=1;
    xtemp = ones(mX,1);
    yrxtemp = yrX;
    L = isnan(X);
    L=(all(L'))';
    xtemp(L)=NaN;
    [xtemp,yrxtemp]=trimnan(xtemp,yrxtemp);
    if ~all([yrX(1)==yrxtemp(1)   yrX(end)==yrxtemp(end)]);
        L= yrX >=yrxtemp(1) & yrX<=yrxtemp(end);
        X=X(L,:);
        yrX=yrxtemp;
        [mX,nX]=size(X);
    end
else;
    [X,yrX]=trimnan(X,yrX);
    [mX,nX]=size(X);
end;
clear xtemp yrxtemp L;

%---- CHECK FOR NEGATIVE CORE INDICES.  FOR ANY SUCH CORE, SCALE THE DEPARTURES FROM THE MEAN
L = X<0;
L1 = any(L) ; % marks cores with any index less than zero
if any(L1);
%     strneg = {'Some core index or residual core index is less than zero for this series',...
%         'Dealing with that by scaling the departures for any such series',...
%         'by the factor mean(x)/(mean(x)-nanmin(x)) '};
%     uiwait(msgbox(strneg,'Message','modal'));
    nflag = sum(L1);
    i1 = find(L1); % xref to cols of X for series to be scaled
    for n = 1:nflag;
        x = X(:,i1(n)); % a core index time series
        xmean=nanmean(x);
        xdep = x-xmean;
        f = xmean / (xmean-nanmin(x)+eps); % the scaling factor for departures
        xnew_dep = f * xdep;
        x1 = xmean+ xnew_dep;
        X(:,i1(n))=x1; % restore the adjusted core index to the tsm
%         if any(x1<0) ; % debug
%             figure(1);
%             plot(yrX,x,yrX,x1)
%             title(['Series' num2str(n)]);
%             grid on;
%             pause
%         end
    end
    if any(any(X<0));
        error('Rescaling still leaves some core index or residual core index negative');
    end
    clear L L1 nflag i1 n x xmean xdep f xnew_dep x1;
 end;

%--- SPECIAL CASE OF JUST ONE CORE INDEX

if nX==1;
    L=isnan(X); % marks which elements of vector X are  NaN
    Result.y=X;
    Result.yry=yrX;
    Result.y1=X;
    Result.yry1 = yrX;
    Result.ylo=X;
    Result.yhi=X;
    Result.ys=zeros(mX,1);
    Result.SSSt = repmat(1,mX,1);
    Result.EPSt = repmat(NaN,mX,1);
    Result.rbar=repmat(NaN,mX,1);
    Result.ntrees=repmat(1,mX,1);
    Result.ncores=repmat(1,mX,1);
    Result.rbs = NaN;
    Result.SSS=1;
    Result.EPS = NaN;
    Result.f = repmat(1,mX,1);
    if any(L);
        Result.ys(L)=NaN
        Result.SSSt(L)=NaN;
        Result.EPSt(L)=NaN;
        Result.rbar(L)=NaN;
        Result.ntrees(L)=0;
        Result.ncores(L)=0;
        Result.f(L)=NaN;
    end;
    return;
end;



%--- HOMOGENIZE VARIANCE OF CORE INDEX-- OPTIONALLY

% 1- compute std dev of each core index
% 2- compute the median of the standard deviations
% 3- adjust each core index to have that standard deviation

if kopt(2)==1;
    xmn = nanmean(X);
    Xmn = repmat(xmn,mX,1);
    X1 = X-Xmn; % core indices as departures from long-term means
    Lmiss = isnan(X1); % keep track of NaN elements in X

    xsd = nanstd(X);
    sd_median = nanmedian(xsd);
    Xsd = repmat(xsd,mX,1);
    Xsd = X1 ./ Xsd; % core indices as standardized departures
    X2 = Xmn + (Xsd * sd_median);  % core indices with mean restored and standard deviation equal to the
    %   median standard deviation of all core indices.
    
    % After the scaling to homogenize variance, of core indices, some core indices may be negative.  Find the
    % largest negative core index and scale departures of core index from 1.0 such that lowest index for any
    % core is 0.
    if any(any(X2<0)); % at least one core index has a negative value
        xmin=nanmin(nanmin(X2)); % minimum index value for any core -- and it is negative
        d=X2-1.0; % departures of core indices from 1.0
        f  = 1 / abs(xmin-1); % the scaling factor for depatures;
        d1 = d*f; % scaled departures
        Xtemp = 1.0 + d1; % restore 1.0 mean
        X2=Xtemp;
        clear Xtemp d1 f d xmin;
    end;
    X2(Lmiss)=NaN; % restore NaN elements
    X=X2;
    clear xmn X1 Lmiss xsd sd_median Xsd X2;
else;
end;



%--- COMPUTE PRELIMINARY SITE CHRONOLOGY-- BEFORE ANY VARIANCE STABILIZATION

% Store time series of lowest and highest core index, and of standard dev of core index
ylo = (nanmin(X'))'; % lowest core index in each year 
yhi = (nanmax(X'))'; % higest ...
ys = (nanstd(X'))'; %  year-by-year standard deviation of core indices

L = ~isnan(X);
n1a = (sum(L'))'; % cv,  number of series in each year
y1=(nanmean(X'))'; % arithmetic average chronology
Lmiss = isnan(y1); % mark missing values in site chron

% Convert site chronology to mean of exactly 1.0
y1mean = nanmean(y1);
y1= (y1 - y1mean) +1.0;
y1(Lmiss)=NaN; % restore missing values
yry1 = yrX;




%--- COMPUTE MEAN BETWEEN-SERIES CORRELATION

[treenms,treexref,check]=treeid(nms); % assign col vector of integer tree ids (treexref)
mwindow=[]; % not using windowed method
tshift =[]; % not using windowed method

% First get the time-invariant r_bar
kopt_r = [3 1 1]; % In call to rbareff, want time invariant rbar, terse mode, and rbar over rows rather than
% cells of correlation matrix
Result=rbareff(X,yrX,minover,kopt_r,treexref,mwindow,tshift);
rsingle=Result.r;
clear Result;

% Then get the pseudo-time dependent r_bar
kopt_r = [2 1 1]; % In call to rbareff, want pseudo-time-dependent rbar, terse mode, and rbar over rows rather than
% cells of correlation matrix
Result=rbareff(X,yrX,minover,kopt_r,treexref,mwindow,tshift);
% .r is time series of mean between-series correlation
% .ncores is time series of sample size (number of cores)
% .npairsn is time series of number of bivariate correlations averaged to get r 
% .nmedian is median length of overlaps of series used for correlation coefficients
n2= Result.ncores;
n1=Result.ntrees;
r=Result.r;


%--- FOR ANY YEARS WITH JUST ONE SERIES, USE LOWEST MEAN BETWEEN-TREE CORRELATION AS THE MEAN BETWEEN-TREE
%   CORRELATION  FOR PURPOSES OF COMPUTING TIME-VARYING SSS AND EPS

% Obsolete section
% L = n2==1;
% rmin = nanmin(nanmin(r));
% rtemp =r;
% rtemp(L)=rmin;
% r=rtemp;
% clear L rtemp rmin;   



%--- SSS AND EPS, time invariant

[sss,epps]=wigley1(rsingle,nanmax(n1));


%--- SSS AND EPS, time pseudo-dependent

[ssst,eppst]=wigley1(r,n1);


%--- VARIANCE ADJUSTMENT -- OPTIONAL

if kopt(1)==1;

    % Compute effective independent sample size, from eq 4, p. 90 of Osborne and Briffa 1992
    c = (n2 - 1) .* r;
    nprime = n2  ./  (1+c);

    % Convert the site index, before adjustment, to departures from long term mean
    ymean =nanmean(y1);  % long term mean
    d1 = y1 - ymean; % index as departures from mean
    L1 = isnan(y1); % mark missing values

    % Scale the departures, from eq 6, p. 92 (time-dependent scaling)
    f=sqrt(nprime);
    d2 = d1 .* f;
    
    % Re-scale the departures to same variance as chronology before variance-stabilization
    ztemp = (d2 - nanmean(d2))/ nanstd(d2);
    utemp = ztemp * nanstd(y1);
    y = ymean + utemp;

   
else;
    y=y1;
    f=repmat(1.0,mX,1);
end;
yry=yry1;


%--- STORE RESULTS

Result.y=y;
Result.yry=yry;
Result.y1=y1;
Result.yry1 = yry1;
Result.ylo=ylo;
Result.yhi=yhi;
Result.ys=ys;
Result.SSSt = ssst;
Result.EPSt = eppst;
Result.rbar=r;
Result.ntrees=n1;
Result.ncores=n2;
Result.rbs = rsingle;
Result.SSS=sss;
Result.EPS = epps;
Result.f = f;