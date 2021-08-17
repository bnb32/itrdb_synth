function Results=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT,specs_OL)
% dtrendrw: convert core ring widths to indices by detrending mean and (optionally) variance
% Results=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT,specs_OL);
% Last Revised 2006-6-3
%
% Convert core ring widths to indices by detrending mean and (optionally) variance.
% Various options available for the detrending,
% including ratio vs difference, RCS vs conventional.
%
%*** IN
%
% W (mX x nW)r ring widths as time series matrix with NaN for missing values
% yrW (mW x 1)i  year vector for W
% idW {? x 1}s ids of cores in W
% yrsW (nW x 2)i  start and end year of period with any data for series in W
% pdfit (1 x 2)i   OPTIONAL start and end year for fit for all series. If [],
%   full length from first to last year of each series is used for the fit. [] must
%   be [] if RCS used
% specs_DT: structure with fields:
% .kopts (2 x 1)i options
%   kopts(1): ratio or difference detrending
%       ==1 ratio
%       ==2 difference
%   kopts(2): emirical growth curve option
%       ==1 %N spline, where N is sample length and % is specified; requests a spline with freq response R
%           at the specified percentage of sample length (typically want this response around 200% the sample
%           length.  % and R are specified in Splinespecs
%       ==2 %N spline, where N is minimum sample length of all series in W, and % specified
%       ==3 %N spline, where N is specified quantile of sample lengths in W, and % specified
%       ==4  n-yr spline, where n specified
%       ==5 monotonically non-increasing spline
%       ==6 horizontal line at mean ring width
%       ==7 modified negative exponential, straight line with neg slope, horizontal line
%           at mean ring width -- in that order of priority
%       ==8 modified negative exponential or straight line of any slope (+ or - or horizontal)
%       ==9 RCS (for future use)
%       ==10 Banded (for future use)
%   kopts(3) variance detrending
%       ==1 detrend variance (remove trend in absolute departures from mean)
%       ==2 do not detrend variance
% .gcrit (1 x 1)r critical value of growth curve; if growth curve, g, drops below this, the series is searched
%       for the longest consecutive segment g id not below gcrit and that segment of index is all that is used
% .yrpith (nW x 1) or []i  pith year of each series (if RCS or banded), or [] (if not RCS or Banded)
% .Splinespecs (nW x 1)r spline information, or [] (see Notes)
% specs_OL == structure of specifications for flagging and removing outlier ring-width values (Notes)
% specs_RCS == structure of specifications if RCS detrending [future use]; [] if N/A
% specs_AB = structure of specifications if age-banded detrending [future use]; [] is N/A
%
%************ OUT
%
% Results -- structure of results, with fields
%
% .W (mW x nW)r  Ring widths as input, except truncated to cover pdfit if pdfit
%       not []
% .yrW (mW x 1)i year vector for W
% .X (mX x nX)r  core indices corresponding to ring widths W; same size as W, except if pdfits not [] --
%   then .X, yrX, etc would only cover pdfits
% .G (mX x nX)r  growth curves, same length as X, with year vector yrX
% .yrX (mX x 1)i year vector for X
% .yrsX (nX x 2)i first and last years of the index (may differ from yrsW if series truncated because
%       of growth curve going to zero or negative
% .ExpVar(nX x 1)r explained variance statistic, computed as
%       1 -   var(e)/var(w),  where e = w - g;  w is ring width, g is fitted curve
% .form (1 x 1)s   "Ratio" or "Difference"
% .period (1 x 1)s  "Full length of ring-width series", "Specified uniform period ????-????"
% .CurveOpt (1 x 1)i  curve option, equivialent to input setting for kopts(2)
% .CurveFinal (nX x 1)i  type of curve actually used for detrending each series;  coded as follows:
%   ==1 spline
%   ==2 modified negative exponential
%   ==3 horizontal line
%   ==4 negative slope straight line
%   ==5 positive slope straight line
%   ==6 specified RCS
% .ParamFit{nX x 1} [variable dims]r  fit parameters; and element of {} for each core
%       If a spline, (2 x 1)r the wavelength of the 0.95 response, and the spline parameter
%       If a modified NE, (3 x 1)r   k,a,b    in  g=  k + a E^(-bt)
%       If a straight line (2 x 1)r  a,b   in    g = a + bt    (where t is 1 at first year of fit)
%       If a horizontal line (1 x 1)r  a, the mean ring width for the fit period
%       If RCS (future use)
%       If Age-Band (future use)
% .explode_flag (nX x 1)L   flag that index shorter than ring width because growth curve went to zero or below
%       at some year;  index before that year is all that is kept
%       ==0  not a problem this series
%       ==1  Yes, this series was truncated
% .Lout {3 x 1}[]L  logical matrics same size as W and X, marking outliers in W. Lout{1} marks row outliers;
%   Lout{3} col outliers, and Lout{3} values that are outliers in rows and cols
% .TrendFraction (nX x 1)r  decimal fraction of variance of ring width that is modeled as trend. Computed as 
%       1 -(var(e)/var(w)), where e is ringwidth minus fitted curve and w is the ringwidth
% .VarDT  Structure of results for variance detrending ([] if kopts(3)==2)
%   .fscale{nX x 1)[? x 1] = tsv of scaling factor for absolute departures for each series 
%   .Numer (nX x 3)r [minimum_fscale maximum_fscale  pctvariance]
%   .explode {}s   running log of any series with problem growth curves -- approaching zero or going negative
% .explode {}s log of any series with problem growth curves in detrending mean ring ring width;
%
%
%*** NOTES
%
% Outlier quality control. The time series matrix of ring widths is analyzed for unusual values in the context
% of a moving time window over individual series, or in the context of all the available widths for a given
% year.  Values flagged as outliers in BOTH categories are optionally converted to NaNs (ring widths and
% indices) in the output matrices W and X.  Results.Lout contains flags for outliers.  See outlier.m for more
% details. 
%
% Caution on outlier removal.  The detrending curves are fit BEFORE any removal of ring-width outliers, and the
% outliers can affect the position of the fitted trend line.  A good strategy if outliers are to be removed is
% first run in diagnostics mode (specs_OT.kremove==0), look at the plots, and remove any outliers from the rwl
% files. Then re-run the detrending. 
%
%
% Splinespecs: if applicable, these depend on kopts(2) setting
%   ==1,2,or 3 specified percentage (e.g., [67  0.95] with kopts(2)==1 means spline with amplitude of
%       frequency response 0.95 at wavelength 67% of length of each series
%   ==4 wavelength (in years) at which spline response is 0.95
%   ==5 [], no info needed
%
% % Specs_OL:  structure with fields:
%   .myr (1 x 1)i   length of window for computing time-varying means and standard deviations (must be odd, e.g., 51) 
%   .nstdev (1 x 2)r  number of standard deviations from mean defining a row outlier (1) and column outlier(2)
%       See outlier.m for more details. Typically set to [3 3]
%   .kremove (1 x 1)L  1 == convert any values that are both row and col outliers to NaN, ==0 just flag them
%       but do not remove
% RCSspecs (future use)
% Bandspecs (future use)
%
% If RCS or Age-Band standardization used, input to this function must include detailed information
% from preliminary steps.  For example, with RCS the detrending curve as a function of year from
% pith must be input
%
% Required that all series in W have at least 30 years of data
%
% Ratio and difference indices.  Both are computed by this function.  The final index output depends on
% kopts(1)
%
% Exploding or negative ratio index. After the ratio index is initially computed, the growth curve g is
% checked for any zero (infinite index) or negative (negative index) values.  The ratio index is then
% truncated to end with the year preceding the first year of zero or negative growth curve. The difference
% index is computed over that same truncated period.
%
% Scaling of index.  The ratio index is shifted to have a mean of exactly 1.0 over the period of the detrending
% fit (after any truncating because of zero or negative growth curve).  The departures from the index from 1.0 are then
% scaled if necessary so that the index is not negative, and so that the index is zero for any year that the ring width is
% zero.  The difference index is scaled to have the same mean and standard deviation as the ratio index, but
% the difference index is not constrained to be positive, and need not be zero in years that the ring width is
% zero.
%
% Variance detrending.  If used, the same spline as for detrending the mean is used to detrend the absoloute
% departures.  After variance detrending, the index is again shifted and scaled such that the mean is 1.0 and
% no negative values occur.  But after variance detrending the correspondence of zero ring width to zero index
% no longer holds


linexp=fittype(@(a,b,c,d,x) a+b*x+c*exp(-d*x));
negexp=fittype(@(a,b,c,x) a+b*exp(-c*x));


%--- CHECK INPUT

[mW,nW]=size(W);
if length(yrW) ~= mW;
    error('yrW not same row size as W');
end;
if ~all(diff(yrW)==1);
    error('yrW does not increment by 1');
end;
L=~isnan(W);
if ~all(sum(L)>30);
    error('All series in W must have at least 30 years of data');
end;

nX=nW;
mX=mW;
X=repmat(NaN,mX,nX); % to hold tsm of core indices
G=repmat(NaN,mX,nX); % to hold tsm of core indices
yrX=yrW;
clear L;


if ~strcmp(class(idW),'cell');
    error('idW must be cell of strings');
end;
if size(idW,2) ~=1;
    error('idW must be col-cell');
end;
if length(idW) ~= nW;
    error('idW must be same length as col-size of W');
end;


[mtemp,ntemp]=size(pdfit);
if ~((mtemp==1 & ntemp==2)  |   isempty(pdfit));
    error ('pdfit must be [] or 1 x 2');
end;
if ~isempty(pdfit);
    if (pdfit(2)-pdfit(1)+1) <= 30;
        error('A specified pdfit must exceed 30 yr');
    end;
    L= yrW>=pdfit(1) & yrW <=pdfit(2);
    Wtemp = W(L,:);
    L = ~isnan(Wtemp);
    if ~all(sum(L')>30);
        error('All series in W must have at least 30 years of data in specified fit period');
    end;
end;
clear L mtemp ntemp Wtemp;


kopts=specs_DT.kopts;
yrpith=specs_DT.yrpith;
Splinespecs = specs_DT.Splinespecs;
gcrit=specs_DT.gcrit;

if ~any(kopts(1)==[1 2])
    error('kopts(1) must be 1 or 2');
end
if ~any(kopts(2)==[1 2 3 4 5 6 7 8 9 10])
    error('kopts(1) must be interger between 1 and 10');
end
if ~any(kopts(2)==[1 2 3 4 8]);
    error('So far, only kopts(2)=[1 2 3 4 8] have been programmed');
end;
if ~any(kopts(3)==[1 2 ])
    error('kopts(3) must be 1 or 2');
end


if ~any(isempty(yrpith) || (size(yrpith,1)==nW && size(yrpith,2)==1));
    error('yrpith must be either empty or nW x 1');
end;

if ~any(    isempty(Splinespecs)  || (size(Splinespecs,2)==2 && size(Splinespecs,1)==1)      );
    error ('Splinespecs must be 1x2  or empty');
end;
if ~isempty(Splinespecs);
    if any(kopts(2) == [ 1 2 3]);
        if Splinespecs(1)<20 || Splinespecs(1)>500;
            error('Pctg of sample length for target wavelength of spline cannot be specified smaller than 20 or larger than 500');
        end;
        if Splinespecs(2)<0.01 || Splinespecs(2)>0.99;
            error('Amplitude of frequency response of spline at target wavelength must be between 0.01 and 0.99');
        end;
    elseif kopts(2)==4;
        if Splinespecs(1)<30 || Splinespecs(1)>(5*nW)
            error('n-year spline cannot be specified with n<30 yr or n>(5 times rowsize of W)');
        end;
        if ~isnan(Splinespecs(2))
            error('Calling for n-year spline, Splinespecs(2) should be NaN');
        end;
    else;
        warning('Splinespecs should be empty if kopts(2) not between 1 and 4');
    end
end;


if ~isa(specs_OL,'struct');
    error('specs_OL must be a structure');
end;
nmstemp =fieldnames(specs_OL);
Sneed = {'myr','nstdev','kremove'};
if ~all(ismember(nmstemp,Sneed));
    error('Not all required fields in specs_OL');
end;
myr = specs_OL.myr;
nstdev = specs_OL.nstdev;
kremove = specs_OL.kremove;

%--- SET SOME SPECS ALSO STORED IN STRUCTURE RESULTS

if kopts(1)==1;
    Results.form='Ratio';
else;
    Results.form ='Difference';
end;

Results.CurveOpt =kopts(2);


% Fit period
if isempty(pdfit)
    Results.period ='Detrending on full length of each ring-width series';
else;
    Results.period =['Detrending restricted to existing data in period ' num2str(pdfit(1)) '-' num2str(pdfit(2))];
end;


%---- INITIALIZE INFO ON PROBLEMATIC GROWTH CURVES

g_explode = {'Problematic Growth Curves'};

%--- OPTIONALLY TRUNCATE DATA TO SPECIFIED RESTRICTED PERIOD FOR FIT

% Truncations due to any specified pdfit
if ~isempty(pdfit);
    L = yrW>=pdfit(1) && yrW<=pdfit(2);
    W=W(L,:);
    X=X(L,:);
    G=G(L,:);
    [mW,nW]=size(W);
    yrW = yrW(L);
    yrX=yrX(L);
    L1 = yrsW(:,1)<pdfit(1);
    L2 = yrsW(:,2)>pdfit(2);
    yrsW(L1,1)=pdfit(1);
    yrsW(L2,1)=pdfit(2);
    clear L L1 L2;
end;
Results.W=W;
Results.yrW=yrW;
Results.yrsX = yrsW;
lensW = yrsW(:,2)-yrsW(:,1)+1;   %   cv of series lengths (could have imbedded NaNs)

%--- FLAGGING OF OUTLIER RING WIDTHS AND OPTIONAL CONVERSION OF THEIR RING-WIDTH  VALUES TO NAN

Lout=outlier(W,yrW,myr,nstdev);
Results.Lout = Lout;
if kremove ==1;
    Lkill = Lout{3};
    W(Lkill)=NaN;
end; %
 
%--- DETRENDING

Results.explode_flag = logical(zeros(nW,1)); % flag for zero or negative growth curve
Results.ParamFit=cell(nW,1);
Results.CurveFinal=repmat(NaN,nW,1);
Results.TrendFraction=repmat(NaN,nW,1);
if kopts(3)==1; % variance detrending option
    Results.VarDT.Numer= repmat(NaN,nW,3); % to store numeric matrix of information
    Results.VarDT.fscale= cell(nW,1); % to hold time series of ratios for variance stabilization
else;
    Results.DT=[];
end;

% Get spline 50% wavelengths
if kopts(2)==1; % spline with specified freq response  at specified % of series length
    pdspline = lensW  .*  (Splinespecs(1)/100); % cv of wavelengths of 0.95 response of desired splines for each series
elseif kopts(2)==2; % base spline wavelength on minimum sample length
    pdspline  = repmat((min(lensW) * (Splinespecs(1)/100)),nW,1); % again  a cv, but all elements same
elseif kopts(2)==3; % base on median sample length
    pdspline  = repmat((median(lensW) * (Splinespecs(1)/100)),nW,1); % again  a cv, but all elements same
elseif kopts(2)==4; % set at specific wavelength for all series
    pdspline  = repmat(Splinespecs(1),nW,1); % again  a cv, but all elements same
else;
    pdspline=[]; % N/A (not a spline) or cannot be set yet (nonincreasing spline)
end;


for n = 1:nW;  % Loop over series
    % Get the time series of ring widths, w/0 leading and trailing NaNs
    w1 = W(:,n); % pull a rw series
    yrw1 = yrW;
    yrgo = yrsW(n,1); % first year of rw series
    yrsp = yrsW(n,2); % last year of rw series
    LW  = (yrW>=yrgo & yrW<=yrsp); % rows of series in W (after optional truncation depending on pdfit)
    w1 = w1(LW); % ring width series with any leading and trailing data sloughed off
    yrw1 = yrw1(LW);

    L1 = isnan(w1); % mark any NaNs (if any, these are internal)
    w1temp = w1; % working copy of ring width time series
    w1temp(L1)=[]; % remove any NaNs from working copy of rw series
    yrw1temp=yrw1;
    yrw1temp(L1)=[]; % year vector for w1temp

 
    if any(kopts(2)== [1 2 3 4]); % if a spline to be fit, and the desired f-response wavelength known already
        f_resp = Splinespecs(2);
        if kopts(2)==4;
            f_resp = 0.5;
        end
        Results.CurveFinal(n)=1; % spline used
        pdthis = pdspline(n); % desired p-response wavelength of spline
        p = splinep(pdthis,f_resp); % spline parameter
        Results.ParamFit{n}=[pdspline(n) p];


        % Call csaps to fit the spline
       
        g = csaps(yrw1temp,w1temp,p,yrw1) ; % cv of spline growth curve, including values at any internal NaNs
        %   in the original ring widths w1

    elseif kopts(2)==8;
        
	Results.CurveFinal(n)=2;
	g = fit(yrw1temp,w1temp,negexp,'StartPoint',[0,1,1]);
	%g = fit(yrw1temp,w1temp,lingexp,'StartPoint',[0,0,1,1]);
	g = g(yrw1temp);
        

    else;
        error('So far programmed only for curve options with kopts(2) one of [1 2 3 4 8]');
    end; %  if any(kopts(2)== [1 2 3 4]); % if a spline to be fit, and the desired 50%-response wavelength known already

    clear w1temp L1 yrw1temp

    %--- Compute the ratio index

    x1 = w1 ./ g; % ratio index
    yrx1 = yrw1;        
    
    %--- TRUNCATE INDEX BEFORE GROWTH CURVE GOES TO ZERO OR BELOW, IF IT EVER DOES
  
    % What if the growth curve approaches zero or goes negative?  If so, pull the longest segment of index for which
    % the growth curve is not so affected; if a tie in length, use earliest segment.
    % Also want to warn user and set an output flag so user knows why the index time
    % series shorter than the ring-width
    L1 = g <gcrit;
    nL1 = sum(L1);
    if nL1>0; % somewhere the growth curve is zero or negative
        yrneg = yrw1(L1);
        yronbad = yrneg(1);
        yroffbad = yrneg(end);
        str_temp = ['  Series ' idW{n}  ': ' num2str(nL1) ' years between ' num2str(yronbad) ' and ' num2str(yroffbad) ' growth curve below ' num2str(gcrit) ' mm'];
        g_explode{length(g_explode)+1}=str_temp;

        gtemp = g;
        gtemp(L1)=NaN;
        L2= longrun(gtemp);
        
        % Store index before truncation, for diagnostic plotting
        xfull = x1;
        yrxfull = yrx1;
        w1full = w1;
        yrw1full = yrw1;
        gfull = g;
               

        x1=x1(L2); % keep just the index for the longest segment with positive growth curve
        w1=w1(L2);
        yrw1=yrw1(L2);
        g=g(L2);
        yrx1 = yrx1(L2);
        
%         % Debug plot
%         figure(1);
%         subplot(2,1,1);
%         plot(yrw1full,w1full,yrw1full,gfull);
%         ylims = get(gca,'YLim');
%         hline = line([yrw1full(1) yrw1full(end)],[gcrit gcrit]);
%         set(hline,'LineStyle','--','Color',[1 0 0]);
%         legend('Ring width','Growth curve','Critical width');
%         grid on;
%         subplot(2,1,2);
%         hp1=plot(yrxfull,xfull,yrx1,x1,[yrxfull(1) yrxfull(end)],[1 1]);
%         set(hp1(1),'Color',[0 0 1]);
%         set(hp1(2),'Color',[1 0 0]);
%         set(hp1(3),'Color',[0 0 0]);
%         pause;
               
        Results.explode_flag(n)=1;
        Results.yrsX(n,:) = [yrx1(1) yrx1(end)];

    end;
    clear L1 L2 nL1 yrneg yronbad yroffbad str_temp str_temp gtemp

    % Force ratio index to mean of exactly 1.0
    doffset = nanmean(x1)-1.0; 
    x1 = x1-doffset;
    
    % With the shifting to exactly mean 1.0 index, it is posssible that ring width of zero does not correspond to index of zero
    % It is also possible that some index is now negative
    L1=w1==0;
    if any(L1);
        dmax = 1.0 - nanmin(x1); % maximum departure of 1.0 from the ratio index;  want this to equal 1.0 since the max departure
        % should occur with the zero ring width
        d = x1 - 1.0; % departures of index from 1.0
        d1 = d * 1.0/dmax; % scaled departures such that lowest index will be 0
        xtemp = 1.0 + d1;
        x1=xtemp;
        L=w1==0;
        if ~all(x1(L)==0);
            error('If ring width is zero, ratio core index should be zero');
        end
        clear dmax d d1 xtemp L;
    elseif any(x1<0);
        dmax = 1.0 - nanmin(x1); % maximum departure of 1.0 from the ratio index;  want this to equal 1.0 since the max departure
        % should occur with the zero ring width
        d = x1 - 1.0; % departures of index from 1.0
        d1 = d * 1.0/dmax; % scaled departures such that lowest index will be 0
        xtemp = 1.0 + d1;
        x1=xtemp;
        clear dmax d d1 xtemp ;
    end;
    if any(x1<0);
        error('Ratio index should not be less than zero');
    end;
    clear L1;



    %------- Compute difference index

    %--- If needed, compute the difference index and substitute it for the ratio index
    % Will scale difference index to same standard deviation as the ratio index. Difference
    % index covers same period as ratio index (i.e., truncated if necessary to include only that
    % part before growth curve g goes to zero or negative
    x1_sd =nanstd(x1); % stdev of ratio index; need this to restore later
    x2=w1-g; % difference index
    yrx2=yrx1;
    
    
    %--- decimal fraction of ring width in fitted trend line  
    
    ftemp =  1-  (nanvar(x2)/nanvar(w1));
    Results.TrendFraction(n)=ftemp;
    
    % Force difference index to mean 1.0 and std dev same as that of ratio index
    x2_sd = nanstd(x2); % stdev of preliminary diff index
    x2_z = (x2-nanmean(x2))/x2_sd; % z-score diff index
    x2 = 1.0 + (x2_z * x1_sd); % scaled difference index: mean 1.0 and standard dev same as that of ratio index
    clear x1_sd x2_sd x2_z


    %------- Use desired form of index as final index

    if kopts(1)==1;
        x = x1;
        yrx=yrx1;
    else;
        x=x2;
        yrx=yrx2;
    end;
    
    
    %--- Optionally detrend variance 
    
    if kopts(3)==1 && any(kopts(2)==[1 2 3 4]);
        kopt_varDT = 2; % not verbose in call to varstab2
        xorig=x;
        yrxorig=yrx;
        [x,yrx,fscale,varexp]=varstab2(x,yrx,p,kopt_varDT);
        % Shift to mean 1.0
        x = x - nanmean(x)+1.0;
        % With variance detrended, what was a zero index may not be a zero index anymore. This is unavoidable,
        % but will want to avoid negative index. Check for any negative index and if found, scale departures of
        % index so lowest value is zero.  
        L = x<0;
        if any(L);
            dmax = 1.0 - nanmin(x); % maximum departure of 1.0 from the ratio index;  want this to equal 1.0 since the max departure
            % should occur with the zero ring width
            d = x - 1.0; % departures of index from 1.0
            d1 = d * 1.0/dmax; % scaled departures such that lowest index will be 0
            xtemp = 1.0 + d1;
            x=xtemp;
            L=w1==0;
            clear dmax d d1 xtemp L;
        end
        clear L;
        Results.VarDT.fscale{n}=fscale;
        Results.VarDT.Numer(n,[1 2 3 ]) = [min(fscale) max(fscale) varexp];
        debugDV=0;

        if debugDV==1;
            plot(yrxorig,xorig,yrx,x);
            grid on;
            legend('Without variance-detrending','With variance-detrending');
            xlabel('year');
        else
        end
    else;
        Results.VarDT=[];
    end;

    %--- SUBSTITUTE COMPUTED INDEX INTO TIME SERIES MATRIX X

    islot = yrx - yrX(1) +1 ; % cv row index
    X(islot,n)=x;
    G(islot,n)=g;
    
    
end; % for n = 1:mW;  % Loop over series


%--- STORE IN RESULTS

Results.W=W;
Results.X=X;
Results.G=G;
Results.yrX = yrX;
Results.explode=g_explode;


