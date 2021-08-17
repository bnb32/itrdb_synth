function [y,yrx,fscale,varexp]=varstab2(x,yrx,p1,kopt)
% varstab2: stabilize variance by spline with specified spline parameter 
% [y,yry,fscale]=varstab1(x,yrx,p1,kopt);
% Last revised 2006-1-3
%
% Stabilize variance by spline with specified spline parameter.
% 
% 
%*** INPUT
%
% x (mx x 1)r time series to be detrended in variance
% yrx (mx x 1)i year vector for x
% p1 (1 x 1)r, spline parameter (for csaps)
% kopt (1 x 1)i  options
%   kopt(1)--- for verbose response with diagnostic plots
%       ==1  Yes for verbose
%       ==2  No: skip the plots
%
%
%*** OUTPUT
%
% y (my x 1)r  series x after detrending variance
% yry (my x 1)i  year vector for y
% fscale (my x 1)r  scaling factor for adjustment of absolute deviations
% varexp (1 x 1)r  percentage of variance of departures explained by fitted trend line
%
%*** UW FUNCTIONS CALLED
%
%
%*** TOOLBOXES NEEDED
%
% spline
%
%
%*** REFERENCES
%
% Meko D. M., Cook E. R., Stahle D. W., Stockton C. W. and Hughes M. K. (1993) Spatial patterns of tree-growth 
% anomalies in the United States and southeastern Canada. J. of Climate 6, 1773-1786. Method described in
% appendix identical except that least-squares straight line used to describe trend.
%
% Cook E. R. and Peters K. (1981) The smoothing spline: a new approach to standardizing forest interior 
% tree-ring width series for dendroclimatic studies. Tree-Ring Bulletin 41, 45-53. Describes relationship 
% between spline parameter a wavelength of 50% response.
%
%*** NOTES
%
% Internal NaNs allowed in the input time series x

% ---- CHECK TIME SERIES INPUT

[mx,nx]=size(x);
if nx~=1; 
    error('x  must be cv');
end;
if mx~=length(yrx);
    error('yrx and row size of x different');
end;
y=repmat(NaN,mx,1); % allocate to store variance-stabilized series


%--- CHECK OPTIONS

if size(kopt,1)~=1 & size(kopt,2)~=1;
    error('kopt should be 1 x 1');
end;
if ~any(kopt(1)==[1 2]);
    error('kopt(1) must be 1 or 2');
end;


%--- TRIM LEADING AND TRAILING NANS

[x1,yrx1]=trimnan(x,yrx);
Lgood = yrx>=yrx1(1) & yrx<=yrx1(end); % pointer to put adjusted series back in original slots of x
if ~all(diff(yrx1)==1);
    error('Year does not increment by 1');
end;
mx1=length(x1);


%-- CHECK  SPLINE PARAMETER p1

if ~isscalar(p1);
else
    if p1<=0 | p1>=1.0;
        error('spline parameter p1 must be between 0 and 1');
    end;
end; 


%-- COMPUTE ABSOLUTE DEVIATIONS

xdev = x1 - nanmean(x1); % deviations from mean
Lnegdev = xdev<0;% pointer to neg deviations
adev = abs(xdev); % absolute deviations




%-- SPLINE FIT TO ABSOLUTE DEVIATIONS
L1 = ~isnan(adev);
xtemp = adev(L1);
yrxtemp =yrx1(L1);
g=csaps(yrxtemp,xtemp,p1,yrx1);
clear xtemp yrxtemp L1;


% REMOVE TREND IN ABS DEVS
fscale = nanmean(adev)./g; % scaling factor is ratio of long-term mean absolute deviation to the trend-line value in each year
bdev=adev .* fscale;  % scaled absolute deviation
bdev(Lnegdev)=-1.0*bdev(Lnegdev); % restore sign to negative departures
xnew = nanmean(x1)+bdev;   % add mean back to departure to get detrended series


% COMPUTE EXPLAINED VARIANCE 
dtemp=adev-g; % difference of absolute departures and trend line
varexp = 100*(1  - (var(dtemp))/var(adev));



% PUT DATA BACK IN ORIGINAL-SIZE VECTOR
y(Lgood)=xnew;

% OPTIONAL

if kopt(1)==1;
    
    figure(1);
    [cL,cB,cW,cH]=figsize(.8,.6);
    set(gcf,'Position',[cL cB cW cH]);
    subplot(2,1,1);
    plot(yrx1,adev,[yrx1(1) yrx1(end)],[mean(adev) mean(adev)],yrx1,g);
    legend('absolute deviations','mean',['spline (' num2str(p1) ' curve']);
    subplot(2,1,2);
    plot(yrx1,x1,yrx1,xnew,[yrx1(1) yrx1(end)],[mean(x1) mean(x1)]);
    legend('Original','Variance Stabilized','Mean');
    grid;
    zoom xon;
    
end;




