function Lout=outlier(X,yrX,myr,nstdev)
% outlier: Flag possible outlier values in a time series vector or matrix
% Lout=outlier(X,yrX,myr,nstdev)
% Last Revised 2006-6-29
%
% Flag possible outlier values in a time series vector or matrix
%
%*** IN
%
% X (mX x nX)r  time series matrix of mX observations, nX series (Notes)
% yrX (mX x 1)i year vector
% myr (1 x 1)i  number of years in the time window for outlier check (Notes--Method); must be odd,
%  and less than the length of any time series in X
% nstdev (1 x 2)r  number of standard deviations from the mean for marking a 
%   row outlier (1) and a column outlier (2).  Flagged as error if less than 3 
%
%*** OUT
%
% Lout {3 x 1}(mX x nX)L   matrices same size as X, with outliers marked as 1, non-outliers as 0, and
%   missing values as NaN. {1} is row-outlier; {2} col outlier; {3} both outlier
%
%*** NOTES
%
% X: may contain missing values coded as NaNs, trailing or ending, or imbedded
%
% Method.   The tsm is checked across rows (years) and variables (columns).
%
% Series are first high-pass filtered with 50-yr spline to avoid identifying as outliers data that are unusually high or low because
% of time trend 
% 
%      Across rows (years) within a column (for an individual time series)
%
%       1 First floor(myr/2) years
%           * myr block beginnning with each year extracted
%           * mean and stdev computed for that part of series
%           * x(i)>= nstdev standard deviations from the mean, year marked as row-outlier
%       2 Last ceil(myr/2) years
%           * same as above, except treating ending block of years
%       3 Other years
%           * myr block centered on key year extracted; data value x(key)
%           * mean and stdev computed for that block
%           * x(key) marked as outlier if necessary following 1 above
% 
%       Across columns (variables)
%           * The time series are converted to zero-mean, unit standard deviation (scaling to offset possible
%               differences in types of variables in the columns (e.g, mix of mountain and desert rainfall records)
%           * row means and standard deviations are computed
%           * any variable >=nstdev standard deviations from row mean is marked as col-outlier for that year
%
%   Any element of the time series matrix that is both a row-outlier and column-outlier is flagged as an
%   outlier

% Hard Code
per = 50; % period for spline detrending
amp = 0.50; % spline will have this freq response at period per
ncrit=5; % a column outlier (outlier within year, across rows) requires at least 5 valid data points in that year

[mX,nX]=size(X);
Lmiss = isnan(X);
La = zeros(mX,nX);
La = logical(La);
Lb = zeros(mX,nX);
Lb = logical(Lb);
Lc = zeros(mX,nX);
Lc = logical(Lc);
Lout=cell(3,1);

%--- CHECK
if myr<=0; myr=1; end

if mod(myr,2)==0;
    error(['myr should be odd, and you have set myr = ' num2str(myr)]);
end
if  any(nstdev<3);
    error(['elements of nstdev should be <3, while you set nstdev = ' num2str(nstdev)]);
end


%--- DETREND SERIES

p = splinep(per,amp);
for n = 1:nX;
    y=X(:,n);
    L=~isnan(y);
    y = y(L);
    x=yrX;
    x=x(L);
    
    s = csaps(x,y,p,x);
%     plot(x,y,x,s);
%     pause;
end;


%--- FLAG ROW OUTLIERS

nends = floor(myr/2); % number of years in the starting and ending blocks
crit_1 = nstdev(1); % number of stnd deviation marking a row-outlier
for n =1:nX; % Loop over series

    x = X(:,n); % get a time series
    yrx = yrX;
    [x,yrx]=trimnan(x,yrx);
    mx = length(x);
    if myr>=mx;
        error('myr must be less than length of series');
    end;
    iX = yrx-yrX(1);
    
      
    %---- Initial Block -- the first  nends years of the series
    
    x1 = x(1:nends);
    yrx1 = yrx(1:nends);
    yrgo_mid=yrx(nends+1);
    irows = yrx1-yrX(1)+1;
    xslab  = x(1:myr); % same slab of years for stats for each test year
    mn1=nanmean(xslab);
    sd1=nanstd(xslab);
      
    L1 = (abs((x1-mn1))>=crit_1*sd1); % col vector marking row-outliers, this block of years, this series
    La(irows,n)=L1;
    clear L1 x1 yrx1 irows xslab mn1 sd1
      
    
    %--- End Block
    
    x1 = x((end-nends+1):end);
    yrx1 = yrx((end-nends+1):end);
    yrsp_mid = yrx(end-nends);
    irows = yrx1-yrX(1)+1;
    xslab  = x((end-myr+1):end); 
    mn1=nanmean(xslab);
    sd1=nanstd(xslab);
    
    L1 = (abs((x1-mn1))>=crit_1*sd1); % col vector marking row-outliers, this block of years, this series
    La(irows,n)=L1;
    clear L1 x1 yrx1 irows xslab mn1 sd1
    
    
    %--- Middle Block (These are centered)
    
    yrxkey=(yrgo_mid:yrsp_mid)'; % cv year vector of values to be tested
    ix = yrxkey-yrx(1)+1; % row index of the key years in x
    xkey = x(ix); % the values to be tested
   
    nmid = length(xkey); % number years in middle block
    iX = yrxkey-yrX(1)+1; % index of test years in X (the target rows)

    delta = floor(myr/2); %half width of myr window (note myr is odd, and half width is width not including central value)
    igo = ix(1)-delta; % row index in x of beginning of slab r xkey(1) testing
    if igo ~=1;
        error('If all was right, igo be 1, as first slab would begin with first data value');
    end;
    
    % Get the segments of time series x for computing stats for testing each value in xkey
    i1 = (igo:(igo+myr-1))'; % cv of row index in x  of first slab for statistics
    I1 = repmat(i1,1,nmid); % dupe cv to matrix
    j1 = (0:(nmid-1));
    J1 = repmat(j1,myr,1);
    I2 = I1+J1; % row index matrix for pulling slabs from x 
    Xslab = x(I2);
    clear i1 I1 j1 J1 I2 ;
        
    % Compute statistics and test xkey for outlier status
    mn1 = (nanmean(Xslab))'; % as col-vector, like xkey
    sd1 = (nanstd(Xslab))';
    L1 = (abs((xkey-mn1))>=crit_1*sd1); % col vector marking row-outliers, this middle block of years, this series
    La(iX,n)=L1;
    clear  L1 mn1 sd1 ;
    
  
        
end; % for n =1:nX; % Loop over series
Lout{1}=La;


%--- FLAG COLUMN OUTLIERS

 crit_2 = nstdev(2); % number of stnd deviation marking a col-outlier


% Convert series to z-scores, col by col
nsize =    (sum((~isnan(X))'))';
mncol = (nanmean(X)); % rv of row means of X
sdcol = (nanstd(X)); % rv of row std devs of X
MNcol = repmat(mncol,mX,1);
SDcol = repmat(sdcol,mX,1);
Z = (X - MNcol) ./ SDcol; % convert to z-scores

% Compute row means of the z-scores
mn1 = (nanmean(Z'))'; % cv of row means of Z
sd1 = (nanstd(Z'))'; % cv of row std devs of Z

% Dupe col vectors of row-means and row-std devs to matrices
Mrow = repmat(mn1,1,nX);
Srow = repmat(sd1,1,nX); 

% Test for outliers
L1 = abs(Z-Mrow) >= crit_2 * Srow;

% If sample size not at least ncrit, do not flag col outliers
L1(nsize<ncrit,:)=0;


Lout{2} = L1;

%--- MARK COMBINED OUTLIERS (BOTH ROW AND COL)

Lout{3} = Lout{1} & Lout{2};

