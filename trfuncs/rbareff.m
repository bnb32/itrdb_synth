function Result=rbareff(X,yrX,minover,kopt,nms,mwindow,tshift)
% rbareff: effective chronology signal for time series matrix
% Result=rbareff(X,yrX,minover,kopt,nms,mwindow,tshift);
% Last revised 2006-6-4
%
% Effective chronology signal for time series matrix. Two or more of the series overlap
% but in general the series may cover different periods.  For tree-ring application, individual
% series may come from different trees, and this tree membership affects the computations.
% The effective chronology signal depends on the sample size (e.g., number of trees) and on the
% interseries correlations.  The signal can optionally be computed in three different ways:
% 1 time-invariant:  correlation matrix based on full overlap of pairs of series
% 2 pseudo time-dependent: correlation matrix based on fulll overlap, but correlations used for a
%   given year restricted to series present
% 3 time-dependent: correlation matrix computed in a sliding time window
%
%*** INPUT
%
% X (mX x nX)r time series matrix (all rows must have at least one non-NaN value);
% yrX (mX x 1)i  time vector (e.g., year) for X
% minover(1 x 1)i  correlation set to NaN if series overlap by fewer than minover years
% nms (nX x 1)i integers assigning tree-membership to series in columns of X (Notes)
% kopt (1 x 2)i  options
%   kopt(1) method
%       ==1 time-dependent: correlations in sliding time window
%       ==2 pseudo time-dependent
%       ==3 time-invariant
%   kopt(2) option for verbose mode (useful for debugging and illustration)
%       ==1 terse mode
%       ==2 verbose mode:  makes graphics; useful in illustrating and
%       debugging
%   kopt(3) correlation use in pseudo time-dependent method (Notes)
%       ==1 mean correlation is average over rows of correlation matrix corresponding to series present
%       ==2 mean correlation is average of correlations for pairs of series present 
% mwindow (1 x 1)i window size. Optional (Notes). Must be odd, cannot exceed 1/3 the
%   row-size of X, but may be empty (Notes).
% tshift (1 x 1)i  shift in windows. Optional (see notes). Must be odd, may be empty (Notes).
%
%
%*** OUTPUT
%
% Result --   structure with fields:
%   .r (mX x 1)r or scalar =  effective chronology signal, or mean between-series correlation in each year
%       (Notes)
%   .yr (mX x 1)i  year vector for ntrees, ncores, npairs, nmedian (Notes)
%   .ncores (mX x 1)i  number of cores in each year
%   .ntrees (mX x 1)i  number of trees in each year ([] is nms [])
%   .npairs (mX x 1)i number of pairs (correlations) mean correlation
%       in each year based on;  scalar if kopt(1)==3
%   .nmedian (mX x 1)r  median length of overlap of series used for mean correlation each year; scalar if kopt(1)==3
%   .eflag (1 x 2) error flags for time-dependent method
%       eflag(1): needed to extrapolate for correlations for last block (one or
%          more ending blocks
%          ==0 no
%          ==1 yes
%   eflag(2): needed to extrapolate for correlations for first block (one or
%          more leading blocks
%          ==0 no
%          ==1 yes
%   eflag=[] -- time invariant method or pseudo time-dependent methods -- blocks not applicable
%
%
%*** REFERENCES
%
% Briffa K. and Jones P. D. (1990) Basic chronology statistics and assessment. In Methods of Dendrochronology,
% Applications in the Environmental Sciences (ed E. R. Cook and L. A. Kairiukstis), pp. 136-152. Kluwer Academic Publishers.
%
% Osborn T. J., Briffa K. R. and Jones P. D. (1997) Adjusting variance for sample-size in tree-ring chronologies
% and other regional mean timeseries. Dendrochronologia 15, 89-99.
%
% Wigley T. M. L., Briffa K. R. and Jones P. D. (1984) On the average value of correlated time series, with
% applications in dendroclimatology and hydrometeorology. Journal of Climate and Applied Meteorology 23, 201-213.
%
%
%*** UW FUNCTIONS CALLED
%
% effcores -- effective number of cores (Briffa and Jones, 1990, p. 142, eq 3.42)
% trimnan
% trailnan
%
%*** TOOLBOXES NEEDED
%
% Stats
%
%*** NOTES
%
% mwindow, tshift.  If kopt(1)==1, analysis is done in sliding time window of width mwindow, shifted tshift
% time units (see below);  If kopt(1)~=2, mwindow and tshift do not apply, are ignored, and may be set to [].
%
% minover:  for any pair of series, correlations are computed, stored, and
% used in the computation of mean between-series correlation only if the
% pair of series overlap by at least minover observations.  For the
% windowed method, the minover coverage must be satisfied within a given
% time window.  For the other methods, the total available overlap of the
% series is used.  If series do not overlap by at least minover years, the corresponding
% bivariate correlation is assigned as NaN
%
% nms:  Computations for tree-ring application require information on which cores (cols of X) belong to
% different trees. This information is passed as input argument nms.
% Best explained by example. Say X has 5 column.  Then nms = [1 1 2 3 3] says columns 1 and 2 are from tree 1,
% column 3 from tree 2, and columns 4 and 5 from tree 3.  For general time series applications,  number the
% columns sequentially (e.g., nms = [1 2 3 4 5]), which is equivalent to saying each series is from a different
% tree. You can set nms manually, or if using tree-ring series following a conventional naming system,
% use function treeid, which  returns the desired nms as output argument treexref.
%
% .r:  The effective chronology signal, a time series or scalar, depending on kopt(1).  r is computed as in
% Briffa and Jones (1990, p. 143, eq 3.43), using the effective number of cores given by eq 3.42 in the same
% reference.  For general time series applications and tree-ring applications with one core per tree, r is
% equivalent to the mean between-series correlation.  If using time-dependent or pseudo time-dependent method,
% r is a time series the same row-size as input X.  For more details on time-dependent r, see 
% Osborn et al. (1997).  If using the time-invariant method, r is a scalar.
% the time-invariant method (kopt(1)==3).
%
% .yr:  if time-invariant method (kopt(1)==3), several of the output arguments (r, npairs, nmedian) are scalars
% rather than time series, and .yr would apply just to .trees and .cores.
%
% kopt(3): In the psuedo time-dependent method (kopt(1)==2), there are two possible ways to include correlations
% for series present in year t in the mean between-series correlation for year t.  First, by kopt(3)==1, is 
% to accept any correlation of a present series with any other series, even if the other series is not present
% in year t.  That method is equivalent to pulling rows of the correlation matrix corresponding to series present.
% Second, by kopt(3)==2 is to accept only correlations between pairs of series present in year t.  That method is
% equivalent to pulling cells of the correlation matrix at intersections of rows and columns for series present in year t.
% If using the time-dependent method (kopt(1)==1), the setting for kopt(3) is ignored and the computations
% proceed as if kopt(3)==2. 
%
% STEPS FOR TIME-DEPENDENT METHOD  (windowed, kopt(1)==1, or Method 1)
%
% 1 - Effective number of cores per tree computed by Wigley and Jones (1990, eqn 3.42).  The equation gives the
%       effective number of cores per tree for each year of the time series.  
% 2- Tsm subdivided into blocks of data of length mwindow shifted by tshift
%       observations.  Blocks have mwindow observations and are aligned so that ending
%       year of last block is last year with data for any of the time series.
%       If first year of the first block does not happen to coincide with the first
%       year with data in any series, an additional leading block, left justified to
%       the first year is formed.
% 3 - Mean between-tree correlation, within-tree correlation, and effective chronology signal
%       computed for each block and assigned to center year of block. Also computed are the number
%       of correlations averaged to get the signal, and the median sample size for correlation
%       computation. No effective chronology signal is computed for a block unless at least two
%       cores from DIFFERENT trees in the block have adequate time coverage for a correlation coefficient.
%       Effective number of cores per tree for the block is set to the median of the annual
%       values of effective number of cores per tree of years in the block.
% 4 - Block-centered values of effective chronology signal, number of pairs, and median sample size
%       linearly interpolated to get estimates at each year of time series.  Before the interpolation,
%       the values for the first and last blocks are assigned to the first and last years of the time series.
% 5 - Interpolated number of pairs set to 0, and median sample size for correlations to NaN for all years with
%       fewer than two trees. This to help user recognize the interpolated effective chronology signal is
%       speculative then.
%
%
%
%
% STEPS FOR PSEUDO-TIME-DEPENDENT METHOD (Method 2)
%
% 1 - Interseries correlation matrix computed using all years of overlap for each pair of series. Separate correlation
%       matrices computed for between-tree and within-tree series.  Correlation set to to NaN if overlap in the
%       two series less than minover years
% 2 - Loop over years of time series matrix X
%   1 identify series (cols of X) present
%   2 compute effective number of cores per tree by eqn 3.42, Briffa and Jones (1990)
%   3 obtain relevant between-tree and within tree correlations (how depends on setting for kopt(3)
%       kopt(3)==1:  use all non-redundant correlations contributed to by any of the series present.  These are
%           the correlations in rows of the correlation matrix corresponding to the series present. If the rows
%           are all-NaN, all non-redundant correlations for the full correlation matrix are used. 
%       kopt(3)==2:  use the nonredundant correlations between pairs of series present.  These are the correlations
%           in cells of the correlation matrix at intersections of rows and columns corresponding to series
%           present.  If those cells are all NaN, all non-redundant correlations for the full correlation matrix are used. 
%   4 compute the average between-tree and within-tree correlation by averaging over the sets of correlations identified
%           above
%   5 compute effective chronology signal by eqn 3.43, Briffa and Jones (1990)
%   6 record the number of correlations contributing to the effective chronology signal.  This is the total of
%           the number of between-tree and within-trees correlations used in computing the mean between-tree
%           and mean within-tree correlation
%   7 record the median sample length used for correlations.  If both between-tree and within-tree correlations
%           used, this is a weigted average of the median sample length of between-tree and within-tree
%           correlations; the weights are proportional to the number of correlations. 
%
%
% STEPS FOR TIME-INVARIANT METHOD (Method 3)
%
% 1 Interseries correlation matrix computed using all years of overlap for each pair of series. Separate correlation
%       matrices computed for between-tree and within-tree series.  Correlation set to to NaN if overlap in the
%       two series less than minover years
% 2 Compute effective number of cores per tree by eqn 3.42, Briffa and Jones (1990)
% 3 Compute average between-tree correlation and average within-tree correlation.  If no valid within-tree correlation
%        (no multiple cores per tree, or none with sufficient overlap), the average within-tree correlation is
%        NaN, and the effective chronology signal reduces to the mean between-tree correlation
% 4 Compute effective chronology signal by eqn 3.43, Briffa and Jones (1990)
% 5 record the number of correlations contributing to the effective chronology signal.  This is the total of
%           the number of between-tree and within-trees correlations used in computing the mean between-tree
%           and mean within-tree correlation
% 6 Record the median sample length used for correlations.  If both between-tree and within-tree correlations
%           used, this is a weigted average of the median sample length of between-tree and within-tree
%           correlations; the weights are proportional to the number of correlations. 
%


%--- CHECK INPUT

% X and yrX
[mX,nX]=size(X);
ntemp=length(yrX);
if mX~=ntemp;
    error('yrX and X not same row size');
end;
yr=yrX;
% % Next section commented out 2006-3-20.  Allow 
%  if any(all(isnan(X')));
%      error('No row of X may be all-NaN');
%  end;


%--- (TRIM TSM TO GET RID OF LEADING AND TRAILING ALL-NAN YEARS)

L1=(all(((isnan(X))')))'; % cv with 1s marking rows in X that are all NaN
if any(L1);
    xtemp=rand(mX,1); % dummy cv, same length as X row size
    xtemp(L1)=NaN;  % row of xtemp is NaN if all series that year NaN in X
    [xtemp,yrtemp]=trimnan(xtemp,yrX); % trim leading and trailing NaNs from xtemp; yrtemp is time vector
    %   for period except the leading and trailing NaNs
    L1=yrX>=yrtemp(1) & yrX<=yrtemp(end);
    X=X(L1,:);
    yrX=yrX(L1);
    [mX,nX]=size(X);
end;
clear L1 xtemp yrtemp;
% (end of trimming)

if kopt(2)==2;
    figure(1);
    plot(yrX,X);
    xlabel('time');
    title(['Plot of time series matrix;, number of series = ' num2str(nX)]);
end



[mtemp,ntemp]=size(kopt);
if ~(mtemp==1 & ntemp==3) | kopt(1)>3 | ~(kopt(2)==1 | kopt(2)==2) | ~(kopt(3)==1 | kopt(3)==2);
    error('kopt must be 1 x 3, kopt(1) must be scalar, and not exceed 3;  kopt(2) and kopt(3)  must be 1 or 2');
end;

if kopt(1)==1;
    [mtemp,ntemp]=size(minover);
    if ~(mtemp==1 & ntemp==1) | minover>yrX;
        error('minover must be scalar, and not greater than row-size of X');
    end;

    [mtemp,ntemp]=size(mwindow);
    if ~(mtemp==1 & ntemp==1) | mwindow>floor(mX/3);
        error('mwindow must be scalar, and not exceed 1/3 row-size of X');
    end;
    kwindow='Yes';
    kmethod=1;
    if mwindow<minover;
        error(['mwindow cannot be smaller than minover']);
    end;
    if tshift>mwindow;
        error('tshift cannot exceed mwindow');
    end;
    if mod(tshift,2)==0 | mod(mwindow,2)==0;
        error('tshift and mwindow must be odd');
    end;
    kopt(3)=2; % for the time-dependent method, require that any correlation used is for pairs of series 
    % present in the block (just having one of the series present will not do). Note that I could remove 
    % this one statement if ever want to have the option to test sensitivity of the windowed method
    % to this restriction
else;
    kwindow='No';
    if kopt(1)==2;
        kmethod=2;
    elseif kopt(1)==3;
        kmethod=3;
    end;
end;


%----- COMPUTE NUMBER OF CORES, TREES AND EFFECTIVE NUMBER OF CORES PER TREE EACH YEAR

[ntrees,ncores,Ceff]=effcores(X,yrX,nms); % implements eq 3.42, Briffa and Jones, 1990

if kopt(2)==2; % debug
    figure(2);
    plot(yrX,ncores,yrX,ntrees);
    legend('cores','trees');
    ylabel('Number of series');
    xlabel('Time');
    grid on;
    title
end
Result.yr = yrX;
Result.ncores=ncores;
Result.ntrees=ntrees;



%--- GET THE TIME SERIES OF CORRELATIONS IN DIFFERENT WAYS

switch kwindow;
    case 'Yes'; % using the time-windowed method (method 1`)

        % ---- MAKE INDEX FOR BLOCK OF ROWS IN X

        btemp=fliplr([mX:-tshift:mwindow]); % index to ending observations of blocks
        ispall=btemp';
        igoall=ispall-mwindow+1;
        if igoall(1)~=1;
            igo1=1;
            isp1=mwindow;
            igoall=[igo1; igoall];
            ispall=[isp1; ispall];
        end;


        % ALLOCATE

        r = repmat(NaN,length(btemp),1); % allocate for correlations for the blocks
        tthis=r;  % allocate for time vecotr for r
        n2 = r; % allocate for number of pairs for each block
        n3=r; % allocate for median length of overlaps for correlations in block

        % LOOP OVER BLOCKS

        i1=0; % initialize counter for blocks
        for n = 1:length(r); % loop over blocks
            i1 = i1+1;
            isp =ispall(n); % ending row of X for this block;
            igo = igoall(n); % start row ...
            tthis(i1) = yrX(1) + igo +  (mwindow-1)/2 -1 ; % time point for middle of block

            X1 = X(igo:isp,:); % subset of data
            yrX1 = yrX(igo:isp);

            % Compute effective number of cores per tree for the block;
            [atemp,btemp,c1]=effcores(X1,yrX1,nms); % c1 is a time series
            c1 =nanmedian(c1); % use the median as representative for the block

            L1=~isnan(X1); % 1 for elements of X1 not NaN
            j1 = (sum(L1)>=minover); % rv: 1 if series has at least minover good data values in block

            % Require that at least two series from DIFFERENT trees have adequate length data in the block.  If
            % not, cannot compute a between-tree correlaton for the block, and will not try to compute an
            % efffective chronology signal

            % Find out whether have enough series
            j2 =length(unique(nms(j1))); % scalar number of series with at least minover years of good data from different trees in block
            if sum(j2)<2;  % % if fewer than 2 series in this time block have required minover time coverage
                r(i1)=NaN;
                n2(i1)=0;
                n3(i1)=NaN;
                continue;
            end;


            % Pull subset of series with at at least minover years of data; no
            % guarantee same minover for each
            X2 = X1(:,j1);
            [mX2,nX2]=size(X2);
            jtree= nms(j1); % cv of tree numbers for series in X2

            % Allocate for correlations between pairs, number of series, and
            %  median overlap of pairs for this time block
            rtemp_b = repmat(NaN,(nX2*(nX2-1))/2,1); % to hold between-tree correlations, at most n*m/2 of these, where m=n-1 and n is number of series
            rtemp_w = repmat(NaN,(nX2*(nX2-1))/2,1); % to hold within-tree correlations, at most n*m/2 of these, where m=n-1 and n is number of series
            n2temp_b=0; % counter for number of between-tree pairs with at least minover overlap
            n2temp_w=0; % counter for number of between-tree pairs with at least minover overlap
            n3temp_b=rtemp_b; % counter for median overlap, between-tree series
            n3temp_w=rtemp_w; % counter for median overlap, between-tree series

            i3_b=0; % between-tree counter
            i3_w=0; % within-tree counter


            % Compute correlations between pairs
            for j3=1:(nX2-1); % outer loop over series
                for j4=(j3+1):nX2; % inner loop over series
                    jtreex=jtree(j3); % tree number for x series
                    jtreey = jtree(j4); % tree number for y series
                    x2=X2(:,j3); % first series of pair
                    y2=X2(:,j4); % second series of pair
                    L2 = ~isnan(x2) & ~isnan(y2);
                    ngood = sum(L2); % number of in-common non-NaN values in the two series
                    if jtreex ~= jtreey; % a between-tree pair
                        i3_b=i3_b+1; % counter specifying target row of rtemp_b, n3temp_b
                        if ngood<minover; % not enough overlap
                            rtemp_b(i3_b)=NaN;
                            n3temp_b(i3_b)=NaN;
                        else;
                            rthis=corrcoef(x2(L2),y2(L2));
                            rthis=rthis(2); % correlation coefficient
                            rtemp_b(i3_b)=rthis;
                            n3temp_b(i3_b)=ngood;
                        end;
                    else; % a within-tree pair
                        i3_w=i3_w+1; % counter specifying target row of rtemp_w, n3temp_w
                        if ngood<minover; % not enough overlap
                            rtemp_w(i3_w)=NaN;
                            n3temp_w(i3_w)=NaN;
                        else;
                            rthis=corrcoef(x2(L2),y2(L2));
                            rthis=rthis(2); % correlation coefficient
                            rtemp_w(i3_w)=rthis;
                            n3temp_w(i3_w)=ngood;
                        end;

                    end ;%  a between-tree pair
                end; % inner loop over series
            end; % outer loop over series

            if all(isnan(rtemp_b));   % If all between-tree correlations for this block are NaN, use NaN for the block
                r(i1)=NaN;
                n2(i1)=0;
                n3(i1)=NaN;
                continue; % go to  next block
            else; % at least one valid between-tree correlation for this block
                % Cull the non-NaN correlations, etc
                rtemp_b = rtemp_b(~isnan(rtemp_b)); % between-tree correlations
                n2temp_b=length(rtemp_b); % number of between-tree correlations
                n3temp_b = n3temp_b(~isnan(n3temp_b)); %  sample sizes for those
                if all(isnan(rtemp_w)) | c1==1; % if no valid within-tree correlations in block, or if effectively 1 core per tree
                    r(i1)=mean(rtemp_b);
                    n3(i1)=median(n3temp_b); % median sample size for correlations
                    n2(i1)=n2temp_b; % number of correlations r(i1) based on
                else; % at least one valid within-tree correlation, and effective number of cores per tree not 1
                    rtemp_w = rtemp_w(~isnan(rtemp_w)); % within-tree correlations
                    n2temp_w=length(rtemp_w); % number of within-tree correlations
                    n3temp_w = n3temp_w(~isnan(n3temp_w)); %  sample sizes for those

                    % Effective chron signal, from Briffa and Jones eq 3.43
                    reff=subfun02(mean(rtemp_b),mean(rtemp_w),c1); % effective chronology signal
                    r(i1)=reff;

                    % Number of correlations contributing to the comptutation of effective chron signal
                    n2(i1)  = n2temp_b + n2temp_w;

                    % Typical sample size for correlations (weighted medians for between-tree and within-tree)
                    n3(i1) = (n2temp_b * median(n3temp_b) + n2temp_w * median(n3temp_w)) / (n2temp_b + n2temp_w);

                end; % if not valid within-tree correlations in block

            end;

        end; % n = 1:tshift: mwindow;

        % Flag if last block does not have enough data for valid
        % correlation
        if isnan(r(end));
            eflag(1)=1;
        else;
            eflag(1)=0;
        end;
        if isnan(r(1));
            eflag(2)=1;
        else;
            eflag(2)=0;
        end;
        
        % SET EFFECTIVE CHRONOLOGY SIGNAL OF FIRST AND LAST YEARS OF THE TIME SERIES  MATRIX TO THE FIRST AND
        % LAST BLOCK SIGNALS
        rstart =r(1); 
        rend = r(end);
        if yrX(1) ~= tthis(1);
            tthis = [yrX(1); tthis];
            r = [rstart; r];
            n2 = [n2(1); n2];
            n3 = [n3(1); n3];
        end
        if yrX(end) ~= tthis(end);
            tthis = [ tthis; yrX(end)];
            r = [r; rend];
            n2 = [n2; n2(end)];
            n3 = [n3; n3(end)];
        end;
        

        %--- LINEARLY INTERPOLATE THE TIME SERIES OF MEAN BETWEEN-SERIES CORRELATION,
        %  SAMPLE SIZE, ETC,TO EACH YEAR OF TIME SERIES MATRIX
        L= ~isnan(r);
        r=interp1q(tthis(L),r(L),yrX); % correlation coeff
        n3=round(interp1q(tthis(L),n3(L),yrX));
        n2=round(interp1q(tthis(L),n2(L),yrX));

        %---- FOR ANY YEARS WITH LESS THAN ONE TREE, REPLACE n2 with zero, n3 with NaN
        n2(ntrees<2) =0;
        n3(ntrees<2) = NaN;
        r(ntrees<2)=NaN;

        %--- STORE RESULTS

        Result.r = r;
        Result.ncores = ncores;
        Result.ntrees=ntrees;
        Result.npairs = n2;
        Result.nmedian=n3;
        Result.eflag=eflag;


        %   .r (mX x 1)r or scalar =  effective chronology signal, or mean between-series correlation in each year
        %       (Notes)
        %   .yr (mX x 1)i, or []  year vector for r (Notes)
        %   .ncores (mX x 1)i  number of cores in each year
        %   .ntrees (mX x 1)i  number of trees in each year ([] is nms [])
        %   .npairs (mX x 1)i number of pairs (correlations) mean correlation
        %       in each year based on;  scalar if kopt(1)==3
        %   .nmedian (mX x 1)r  median length of overlap of series used for mean correlation each year; scalar if kopt(1)==3
        %   .eflag (1 x 2) error flags
        %       eflag(1): needed to extrapolate for correlations for last block (one or
        %          more ending blocks
        %          ==0 no
        %          ==1 yes
        %   eflag(2): needed to extrapolate for correlations for first block (one or
        %          more leading blocks
        %          ==0 no
        %          ==1 yes
        %   eflag=[] -- time invariant method -- blocks not applicable
        %
    otherwise; % not using sliding window;  either kopt==2 or kopt==3;


        %------ COMPUTE BETWEEN-SERIES CORRELATIONS, USING AVAILABLE OVERLAP


        R=repmat(NaN,nX,nX); % allocate for correlation matrix
        Tdiff = logical(ones(nX)); % initialize as if all series were from different trees
        NR  = repmat(NaN,nX,nX); % allocate for sample size (# obs) for each correlation in R

        j=1:nX;
        for n =1:nX; % outer loop over series
            %         if n==2; % debug
            %             disp('here');
            %         end;

            tree_outer = nms(n);
            x1= X(:,n);
            j1 = j ~= n;
            i1=find(j1);
            ninner =sum(j1);

            for m = 1:ninner; % inner loop over series
                mthis =i1(m);
                tree_inner = nms(mthis);
                y1=X(:,mthis);

                % --- Get overlap of x1 and y1

                L1= ~isnan(x1) & ~isnan(y1);
                nover = sum(L1);
                if nover<minover; % insufficient overlap for acceptable correlation
                    continue;
                else;
                    x2=x1(L1);
                    y2=y1(L1);
%                     if all(isnan(x2)) | all(isnan(y2));
%                         disp('here')
%                     end
                    rthis=corrcoef(x2,y2);

                    if tree_outer==tree_inner;
                        Tdiff(n,mthis)=0; % mark this correlation as within tree
                    end;
                    R(n,mthis)=rthis(2); % store correlation coefficient
                    NR(n,mthis)=length(x2); % store sample size for the correlation
                end;
            end; % inner loop over series
        end; %  outer loop over series

        % Make copies of correlation matrix applicable to using between-tree and within-tree correlations only
        RB=R; % correlations, between-tree
        RW=R; % within-tree
        NRB=NR; % median sample size for correlations
        NRW=NR;
        RB(~Tdiff)=NaN;
        RW(Tdiff)=NaN;
        NRB(~Tdiff)=NaN;
        NRW(Tdiff)=NaN;



        if kmethod==2; % time varying, in a way
            r=repmat(NaN,mX,1); % to hold the time series of effective chronology signal
            n2=repmat(NaN,mX,1); % to hold the number of correlations that mean based on
            n3=repmat(NaN,mX,1); % to hold median sample size of time series pairs for correlations
            Lthis = logical(zeros(mX,nX));
            for n = 1:mX; % loop over years
                ceff=Ceff(n);
                xthis = X(n,:);
                L1 = ~isnan(xthis);
                Lthis(n,:)=L1;
                nuse=sum(L1); % use all series with valid data this year
                if ceff==1 | all(isnan(RW(:))); % if effective number of cores per tree is 1 or if no within-tree correlations
                    % must rely on between-tree correlation matrix only
                    k=1;
                else
                    k=2;
%                     if n==93;
%                         disp('here')
%                     end
                end;
                [r(n),n2(n),n3(n)]=subfun03(RB,RW,ncores(n),NRB,NRW,ceff,L1,k,kopt(3));
                

                eflag=[];
            end; %  for n = 1:mX; % loop over years
            Result.r=r; % effective chronology signal



        elseif kmethod==3; % time invariant

            ceff=nanmedian(Ceff); % set a global effective number of cores per tree
            % By using median rather than mean, avoid overly high weight on multiple-core sampling
            % usually found with the younger trees.   It is the early part of the record most relevant to loss
            % of signal from decreasing sample size.

            L1 = logical(ones(1,length(nms))); % indicates that all series used in computation

            if ceff==1 | all(isnan(RW(:))); % if effective number of cores per tree is 1 or if no within-tree correlations
                % must rely on between-tree correlation matrix only
                k=1;
            else
                k=2;
            end;
            nc = length(nms); % number of cores is total number of cols of X
            [r,n2,n3]=subfun03(RB,RW,nc,NRB,NRW,ceff,L1,k,kopt(3));
            Result.r=r; % effective chronology signal

            eflag=[];
        end;
end; %  switch kwindow


Result.npairs = n2; % number of correlations used in computation of effective chronology signal
Result.nmedian = n3; % median length of overlap used for correlations
Result.eflag=eflag; % error flag relevant only if using windowed method



%---------------------SUBFUNCTIONS

function J1=subfun01(j1)
%
% Builds row and col indices to a reference correlation matrix
% j1 are the col indices (variables) whose elements of matrix (row, col)
% are desired
% J1 are the rows and cols to pull from the correlation matrix
% Redundant row/col combinations are excluded: (e.g., row 4/col 3 and row 3/ col 4 are redundant, only one is
% used)

nsize =length(j1);
for n=1:(nsize-1);
    if n==1;
        i1=repmat(j1(1),(nsize-1),1);
        i2 = (j1(2:nsize))';
    else;
        i1=[i1; repmat(j1(n),(nsize-n),1)];
        i2=[i2; (j1((n+1):nsize))'];
    end;
end;
J1=[i1 i2];

%-----------------------------------

function reff=subfun02(rbt,rwt,ceff)
%
% Compute the effectivive chronology signal from the mean between-tree correlation, mean
% within-tree correlation, and effective number of cores per tree
%
% 2006-6-4 modification.  Found a chronology AZ103.rwl for which r_bt>r_wt, resulting in
% a computed reff>1.0.  This result is illogical.  To avoid reff>1.0, the function now checks
% for rbt>rwt, and, if true, sets reff to the lower of rbt and rwt.
%
% Briffa and Jones, 1990, eqn 3.43
if ceff==1;
    reff=rbt;
elseif rbt>=rwt;
    reff = rwt; 
else
    reff = rbt /  (rwt + (1-rwt)/ceff);
end

%----------------------------------

function [r,nr,nmed]=subfun03(RB,RW,nc,NRB,NRW,ceff,L1,k,kopt3)
% Computation of the mean between-tree and within-tree correlations for use in the
% time-invariant and pseudo time-dependent methods for effective chronolgy signal.
% Computations apply for a single year (if pseudo time-dependent method)
% RB = correlation matrix of between-tree series
% RW = correlation matrix of within-tree series
% nc = number of cores this year
% NRB = matrix of number of observations for correlations in RB 
% NRW = matrix of number of observations for correlations in RW
% ceff = effective number of cores per tree this year
% L1 = logical pointer to series present this year
% k = 1 means either no multiple cores per tree or no non-NaN within-tree correlations; must use RB, NRB for all calculations
% kopt3 = option for using rows vs cells of correlation matrix

rowkey = find(L1);
Rkey = RB(rowkey,:);
NRkey = NRB(rowkey,:);

if k==1; % either no multiple cores per tree or no non-NaN within-tree correlations; must use RB, NRB for all calculations
    if all(isnan(Rkey(:))); % no valid correlations for any this-year series
        r=nanmean(RB(:));
        nr =   sum(~isnan(RB(:)))/2;
        nmed = round(nanmedian(NRB(:)));
    else; % at least one valid correlation for a this-year series, though not necessarily with another this-year series
        % Find out if any valid bivariate correlations for at least one pair of this-year cores
        if nc ~=1;
            J1=subfun01(find(L1));
            irow=J1(:,1);
            icol=J1(:,2);
            ilinear=sub2ind(size(RB),irow,icol); % gives single-col index into correl matrix
        else; 
        end;
        if (nc~=1) & (~all(isnan(RB(ilinear))));
            if kopt3==2;
                r=nanmean(RB(ilinear));
                nr = sum(~isnan(RB(ilinear)));
                nmed = round(nanmedian(NRB(ilinear)));
            else;
                RB1 = RB;
                NRB1 = NRB;
                RB1(ilinear)=NaN;
                NRB1(ilinear)=NaN;
                Rkey = RB1(rowkey,:);
                NRkey = NRB1(rowkey,:);
                r=nanmean(Rkey(:));
                nr = sum(~isnan(Rkey(:)));
                nmed = round(nanmedian(NRkey(:)));
            end;
        else;
            r=nanmean(Rkey(:));
            nr = sum(~isnan(Rkey(:)));
            nmed = round(nanmedian(NRkey(:)));
        end;
    end;
else; % can possibly use within-tree correlations
    
    Rwkey = RW(rowkey,:);
    NRwkey = NRW(rowkey,:);
   
    %--- Between-tree part of computation
    
    if all(isnan(Rkey(:))); % no valid between-tree correlations for any this-year series
        r_bt=nanmean(RB(:));
        nr_bt =   sum(~isnan(RB(:)))/2;
        nmed_bt = nanmedian(NRB(:));
    else; % at least one valid between-tree correlation for a this-year series
        % Find out if any valid bivariate correlations for at least one pair of this-year cores
        if nc~=1;
            J1=subfun01(find(L1));
            irow=J1(:,1);
            icol=J1(:,2);
            ilinear=sub2ind(size(RB),irow,icol); % gives single-col index into correl matrix
        else
        end
        if (nc~=1) & (~all(isnan(RB(ilinear))));
            if kopt3==2;
                r_bt=nanmean(RB(ilinear));
                nr_bt = sum(~isnan(RB(ilinear)));
                nmed_bt = nanmedian(NRB(ilinear));
            else;
                RB1 = RB;
                NRB1 = NRB;
                RB1(ilinear)=NaN;
                NRB1(ilinear)=NaN;
                Rkey = RB1(rowkey,:);
                NRkey = NRB1(rowkey,:);
                r_bt=nanmean(Rkey(:));
                nr_bt = sum(~isnan(Rkey(:)));
                nmed_bt = round(nanmedian(NRkey(:)));
            end
        else;
            r_bt=nanmean(Rkey(:));
            nr_bt = sum(~isnan(Rkey(:)));
            nmed_bt = nanmedian(NRkey(:));
        end;

    end;

    %--- Within-tree part of computation
    
    if all(isnan(Rwkey(:))); % no valid within-tree correlations for any this-year series
        r_wt=nanmean(RW(:));
        nr_wt =   sum(~isnan(RW(:)))/2;
        nmed_wt = nanmedian(NRW(:));
    else; % at least one valid within-tree correlation for a this-year series
        % Find out if any valid bivariate correlations for at least one pair of this-year cores
        if nc~=1;
            J1=subfun01(find(L1));
            irow=J1(:,1);
            icol=J1(:,2);
            ilinear=sub2ind(size(RW),irow,icol); % gives single-col index into correl matrix
        else;
        end;
        if (nc~=1) & (~all(isnan(RW(ilinear))));
            if kopt3==2;
                r_wt=nanmean(RW(ilinear));
                nr_wt = sum(~isnan(RW(ilinear)));
                nmed_wt = nanmedian(NRW(ilinear));
            else
                RW1 = RW;
                NRW1 = NRW;
                RW1(ilinear)=NaN;
                NRW1(ilinear)=NaN;
                Rkey = RW1(rowkey,:);
                NRkey = NRW1(rowkey,:);
                r_wt=nanmean(Rkey(:));
                nr_wt = sum(~isnan(Rkey(:)));
                nmed_wt = round(nanmedian(NRkey(:)));
            end
        else;
            r_wt=nanmean(Rwkey(:));
            nr_wt = sum(~isnan(Rwkey(:)));
            nmed_wt = nanmedian(NRwkey(:));
        end;
    end;

    %--- Effective chronology signal, and associated data
    
    % For the median sample size of correlations, use the weighted median sample sizes of pairs used for
    % between-tree and within-tree mean correlation.  The weights are proportional to the number of
    % correlations used.  Final value rounded
    nmed = round((nmed_bt * nr_bt + nmed_wt * nr_wt)/(nr_bt + nr_wt));

    % For the number of correlatons used for the effective chronology signal, use the sum of the number of
    % between-tree and within-tree correlations
    nr = nr_bt + nr_wt;

    % Compute effective chronology signal
    reff=subfun02(r_bt,r_wt,ceff);
    r=reff;

end
