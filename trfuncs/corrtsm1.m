function Result=corrtsm1(X,yrX,minover,nms,kopt,mwindow,tshift)
%function [r,yrr,n1,n2,n3,ntrees,eflag]=corrtsm1(X,yrX,minover,kopt,mwindow,tshift,nms)
% corrtsm1: mean between-series or between-tree correlation for a time series matrix
% Result=corrtsm1(X,yrX,minover,nms,kopt,mwindow,tshift);
% Last revised 2006-2-23
%
% Mean between-series or between-tree correlation for a time series matrix. Given a time
% series matrix whose series overlap somewhere but may cover different time periods,
% compute the mean correlation between series.  This mean correlation is a scalar
% if time-dependence is not considered and a vector if time-dependence is
% considered.  Two types of time-dependence are optionally handled:
% 1 true time dependence: correlations are computed in a sliding time window
% 2 pseudo time-dependence:  between series correlations computed from full overlap of
%   the time series, but mean between-series correlation for a given year computed using only
%   correlations for series with data in that year
% See notes for more.
%
%
%*** INPUT
%
% X(mX x nX)r time series
% yrX (mX x 1)i  time vector (e.g., year) for X
% minover(1 x 1)i  minimum overlap in a pair of series for acceptable
%   correlation coefficient
% nms (nX x 1)i integers associating each column of X with a different "tree" (Notes)
% kopt (1 x 2)i  options
%   kopt(1) method
%       ==1 correlations in sliding time window
%       ==2 correlation as time for year t from full overlap of pairs of series with data in year t
%       ==3 correlation for year t is mean of correlations for pairs of series,
%           regardless if they have data in year t; correlations by this
%           option not time dependent
%   kopt(2) option for verbose mode
%       ==1 terse mode
%       ==2 verbose mode:  makes graphics; useful in illustrating and
%       debugging
% mwindow (1 x 1)i window size. Optional(see notes). Must be odd, cannot exceed 1/3 the
%   row-size of X, but may be empty (see notes).
% tshift (1 x 1)i  shift in windows. Optional (see notes). Must be odd, may be empty (see notes).
% nms {? x 1}s or []:  for general time series, [] applies.  If tree-ring application and want  
%   between-tree correlation, must tell which series belong to which trees.  In that case, nms
%   holds the site-tree-core ids (Notes)
%
%*** OUTPUT
%
% r (? x 1)r  mean between-series correlation in each year; ? may be less
%   than or equal to mX; less than if fewer than 2 series in early years
% yrr (? x 1)i  year vector for r
% n1 (? x 1)i  number of series with data each year
% n2(? x 1)i or (1 x 1)i interpolated number of pairs mean correlation
%   based on;  scalar if kopt(1)==3
% n3 (? x 1)r  median length of overlap of series correlations computed
%   on; scalar if kopt(1)==3
% ntrees (? x 1)i  number of trees the mean between-tree correlation r is based on
%   [] returned if nms is [].
% eflag (1 x 2) error flags
%   eflag(1): needed to extrapolate for correlations for last block (one or
%       more ending blocks
%       ==0 no
%       ==1 yes
%   eflag(2): needed to extrapolate for correlations for first block (one or
%       more leading blocks
%       ==0 no
%       ==1 yes
%   eflag=[] -- time invariant method -- blocks not applicable
%
%*** REFERENCES -- none
%
%*** UW FUNCTION CALLED
%
% effcores
% trimnan -- trims leading and trailing all-NaN rows from X
% treeid
%
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% mwindow, tshift.  If kopt(1)==1, analysis is done in sliding time window of width mwindow, shifted tshift
% time units (see below);  If kopt(1)~=2, mwindow and tshift do not apply, and should be set to [].
%
% minover:  for any pair of series, correlations are computed, stored, and
% used in the computation of mean between-series correlation only if the
% pair of series overlap by at least minover observations.  For the
% windowed method, the minover coverage must be satisfied within a given
% time window.  For the other methods, the total available overlap of the
% series is used.  If series do not overlap by at least minover years, the corresponding
% bivariate correlation is assigned as NaN
%
% nms:  Because main purpose of function was to compute mean between-tree and between-core correlation
% for time series matrix of tree-ring indices, need to indicate which trees the time series in X belong to. 
% Best explained by example. Say X has 5 column.  Then nms = [1 1 2 3 3] means columns 1 and 2 are from tree 1,
% column 3 from tree 2, and columns 4 and 5 from tree 3.  For general time series applications, can set 
% nms = [], meaning each column represents a separate tree.  The same could be be indicated by numbering the
% columns sequentially (e.g., nms = [1 2 3 4 5].  You can set the values for nms manually before calling
% corrtsm1, or use function treeid, which would return the desired nms as output argument treexref.
%
% STEPS FOR WINDOWED METHOD (true time-dependence, or Method 1)
% 1- Tsm subdivided into blocks of data of length mwindow shifted by tshift
%   observations.  Blocks have mwindow observations and are aligned so that ending
%   year of last block is last year with data for any of the time series.
% 2 - Mean between-tree correlation computed for each block and assigned to end
%   year of block.
% 3 - If first year of the first block does happen to coincide with the first
%   year with data in any series, an additional leading block, left justified to
%   the first year is formed, and the correlations computed for that block
% 4 - Linear interpolation is finally used to assign mean-correlation values to each of the time series
%       matrix.  The correlations for each block are first assigned to the middle years of the block.
%       The last available correlation is assigned to the last year of the time series, the first block-correlation
%       to the first year of the time series.  Correlation is then linearly interpolated for all years in
%       between.
% 5 - Number of pairs (correlations) for each block are counted and stored and similarly
%       interpolated. These are rounded to integers
% 6 - Median number of years of overlap for series within each block
%       similarly counted and stored as integers
% 7 - Special flags eflag used to warn user that needed to either first or
%     last block correlations also needed to be interpolated
% 8 - Special conditions apply for leading or trailing periods with just one time series (or, if nms not [], just one tree).
%   In that case cannot compute correlations.  The first available valid correlation is assigned to the leading years, the last
%   available correlation to the trailing years.  n2 and n3 are set to zero for those years.  
%
% STEPS FOR NON-WINDOWED TIME-VARIABLE-COMPOSITION METHOD (Method 2)
%
% 1 - Interseries correlations computed using all years of overlap for each pair;  if nms not [], only the correlations
%       between series from different trees are used in later computations 
% 2 - Series with data in each year identified
% 3 - Pointer from (2) used to pull relevant correlations
% 4 - Mean correlation computed for each year using correlations for series
%       with data in that year
% 5 - n1 is a time series, same length as r
% 6 - n2 and n3 are time series, also same length as r, but generally different
%       than by method 1: a series can contribute to the correlations as
%       long as it has data in the given year and at some time in the
%       record has at least  minover years overlap with another of the
%       series with data in that year.  The median length of overlap now not constrained by
%       mwindow, so may be larger than by method 1.
%
% STEPS FOR TIME-INVARIANT METHOD (Method 3)
%
% 1 - Interseries correlations computed using all years of overlap for each pair
% 2 - Correlations averaged over pairs (with constraint that a correlation
%       requires minover years overlap in series
% 3 - returned r is a scalar because r not variable with time
% 4 - Similarly, n2 and n3 scalars


%--- CHECK INPUT


[mX,nX]=size(X);
ntemp=length(yrX);
if mX~=ntemp;
    error('yrX and X not same row size');
end;
yrr=yrX;


%---- GET THE BETWEEN-TREE INFORMATION

if isempty(nms);
    treenms=[];
    treexref=[];
    check=[];
else;
    [treenms,treexref,check]=treeid(nms);
end;



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
    yrr=yrX;
end;
clear L1 xtemp yrtemp;
% (end of trimming)

if kopt(2)==2;
    figure(1);
    plot(yrX,X);
    xlabel('time');
    title(['Plot of time series matrix;, number of series = ' num2str(nX)]);
end

[mtemp,ntemp]=size(minover);
if ~(mtemp==1 & ntemp==1) | minover>yrX;
    error('minover must be scalar, and not greater than row-size of X');
end;

[mtemp,ntemp]=size(kopt);
if ~(mtemp==1 & ntemp==2) | kopt(1)>3 | ~(kopt(2)==1 | kopt(2)==2);
    error('kopt(1) must be scalar, and not exceed 3;  kopt(2) must be 1 or 2');
end;

if kopt(1)==1;
    [mtemp,ntemp]=size(mwindow);
    if ~(mtemp==1 & ntemp==1) | mwindow>floor(mX/3);
        error('mwindow must be scalar, and not exceed 1/3 row-size of X');
    end;
    kwindow='Yes';
    kmethod=1;
else;
    if ~(isempty(mwindow) & isempty(tshift));
        error('mwindow and tshift must be [] if kopt~=1');
    else;
    end;
    kwindow='No';
    if kopt(1)==2;
        kmethod=2;
    elseif kopt(1)==3;
        kmethod=3;
    end;
end;


%----- COMPUTE NUMBER OF SERIES IN EACH YEAR

n1 =   (sum(~isnan(X')))';

if kopt(2)==2;
    figure(2);
    plot(yrX,n1);
    ylabel('Number of series');
    xlabel('Time');
    title('Sample size in each year');
end;


%----- OPTIONALLY COMPUTE NUMBER OF TREES AND EFFECTIVE NUMBER OF CORES PER TREE EACH YEAR 

if ~isempty(nms);
    [ntrees,ncores,Ceff]=effcores(X,yrX,nms);
    clear ncores;
else;
    ntrees=[];
    Ceff=[]
end


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

            L1=~isnan(X1); % 1 for elements of X1 not NaN
            j1 = (sum(L1)>=minover); % rv: 1 if series has at least minover good data values in block
            
            % Find out whether have enough series or trees for at least one correlation coefficient
            if ~isempty(nms);
                j2 =length(unique(treexref(j1))); % scalar, n of unique trees in block with at least minover good data
            else;
                j2=sum(j1);  % scalar=number of series in block with at least minover good data values
            end;
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
            
            if ~isempty(nms);
                jtree= treexref(j1); % cv of tree numbers for series in X2
            end;

            % Allocate for correlations between pairs, number of series, and
            %  median overlap of pairs for this time block
            rtemp = repmat(NaN,(nX2*(nX2-1))/2,1); % to hold correlations, at most n*m/2 of these, where m=n-1 and n is number of series
            n2temp=0; % counter for number of pairs with at least minover overlap
            n3temp=rtemp; % counter for median overlap

            i3=0;


            % Compute correlations between pairs, optionally ignoring any within-tree correlations
            for j3=1:(nX2-1);
                for j4=(j3+1):nX2;
                    if ~isempty(nms);
                        jtreex=jtree(j3); % tree number for x series
                        jtreey = jtree(j4); % tree number for y series
                    else;
                        jtreex=1;
                        jtreey=2; % these dummy values guarantee series considered from different trees
                    end;
                    i3=i3+1; % counter specifying target row of rtemp, n3temp
                    x2=X2(:,j3); % first series of pair
                    y2=X2(:,j4); % second series of pair
                    L2 = ~isnan(x2) & ~isnan(y2);
                    ngood = sum(L2); % number of in-common non-NaN values in the two series
                    if ngood<minover;
                        % not enough overlap
                        rtemp(i3)=NaN;
                        n3temp(i3)=NaN;
                    else;
                        rthis=corrcoef(x2(L2),y2(L2));
                        rthis=rthis(2);

                        rtemp(i3)=rthis;
                        n3temp(i3)=ngood;

                        % Correlation not to be used if nms not [] and series from same tree
                        if ~isempty(nms);
                            if jtreex==jtreey;
                                rtemp(i3)=NaN;
                                n3temp(i3)=NaN;
                            end;
                        end
                    end;
                end;
            end;



            % If all correlations for this block are NaN, use NaN for the block
            % summary
            if all(isnan(rtemp));   % If all correlations for this block are NaN, use NaN for the block
                % summary
                r(i1)=NaN;
                n2(i1)=0;
                n3(i1)=NaN;
                continue;
            else
                if any(isnan(rtemp)); % if any of the correlations are NaN, do not use them
                    rtemp=rtemp(~isnan(rtemp));
                    n3temp=n3temp(~isnan(n3temp));
                else;
                end;
            end;

            % Fill vectors
            r(i1)=mean(rtemp);
            n3(i1)=median(n3temp);
            n2(i1)=length(rtemp);

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
        end; % n = 1:tshift: mwindow;


        % ASSIGN COMPUTED CORRELATION, ETC, OF FIRST AND LAST BLOCKS TO THE
        % FIRST AND LAST YEARS OF THE TIME SERIES MATRIX
        tthis=[yrX(1); tthis; yrX(end)];
        r = [r(1) ; r; r(end)];
        n2 = [n2(1) ; n2; n2(end)];
        n3 = [n3(1) ; n3; n3(end)];


        %--- LINEARLY INTERPOLATE THE TIME SERIES OF MEAN BETWEEN-SERIES CORRELATION,
        %  SAMPLE SIZE, ETC,TO EACH YEAR OF TIME SERIES MATRIX
        r=interp1q(tthis,r,yrX); % correlation coeff
        n3=round(interp1q(tthis,n3,yrX));
        n2=round(interp1q(tthis,n2,yrX));
        
        %---- REPLACE ANY REMAINING LEADING NANS WITH FIRST r, n2=0, n3=0
        itemp=find(~isnan(r));
        if itemp(1) ~= 1;
            r(1:itemp(1))=r(itemp(1));
            n2(1:(itemp(1)-1)) = 0; 
            n3(1:(itemp(1)-1)) = 0;
        end
        
        
        %---- REPLACE ANY TRAILING NANS WITH LAST r, n2=0, n3=0
        itemp=find(~isnan(r));
        if itemp(end) ~= length(r);
            r((itemp(end)+1):end) =r(itemp(end));
            n2((itemp(end)+1):end) =0;
            n3((itemp(end)+1):end) =0;
        end
        
        

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
            if ~isempty(nms);
                tree_outer = treexref(n);
            else
                tree_outer=[];
            end;
            x1= X(:,n);
            j1 = j ~= n;
            i1=find(j1);
            ninner =sum(j1);

            for m = 1:ninner; % inner loop over series
                mthis =i1(m);
                if ~isempty(nms);
                    tree_inner = treexref(mthis);
                else
                    tree_inner=[];
                end;
                y1=X(:,mthis);

                % --- Get overlap of x1 and y1

                L1= ~isnan(x1) & ~isnan(y1);
                nover = sum(L1);
                if nover<minover; % insufficient overlap for acceptable correlation
                    continue;
                else;
                    x2=x1(L1);
                    y2=y1(L1);
                    rthis=corrcoef(x2,y2);
                    if isempty(nms);
                    else
                        if tree_outer==tree_inner;
                            Tdiff(n,mthis)=0; % mark this correlation as within tree
                        end;
                    end;
                    R(n,mthis)=rthis(2); % store correlation coefficient
                    NR(n,mthis)=length(x2); % store sample size for the correlation
                end;
            end; % inner loop over series
        end; %  outer loop over series

        if kmethod==2; % time varying, in a way
            
            % Make copies of correlation matrix applicable to using between-tree and within-tree correlations only
            RR=R; % between-tree
            RW=R; % within-tree
            NRR=NR;
            if ~isempty(nms);
                RR(~Tdiff)=NaN;
                RW(Tdiff)=NaN;
                NRR(~Tdiff)=NaN;
            end;
            
            r=repmat(NaN,mX,1); % to hold the time series of mean between-series r
            n2=repmat(NaN,mX,1); % to hold the number of correlations that mean based on 
            n3=repmat(NaN,mX,1); % to hold median sample size of time series pairs for correlations
            Lthis = logical(zeros(mX,nX));
            for n = 1:mX; % loop over years
                ceff=Ceff(n);
                xthis = X(n,:);
                L1 = ~isnan(xthis);
                Lthis(n,:)=L1;
                nuse=sum(L1); % use all series with valid data this year
                
                if nuse==1; % only one series available; use as r the average between-tree correlation of that series with 
                    % all others.  If the series has a NaN as its correlation matrix entry, use the overall
                    % average of the between-tree correlations as r
                    rowkey = find(L1);
                    Rkey = RR(rowkey,:);
                    if all(isnan(Rkey)); 
                        r(n)=nanmean(RR(:));
                    else;
                        r(n)=nanmean(Rkey);
                    end;
                    
                elseif nuse>1 & ntrees(n)==1; % just one tree, and multiple cores from it; use eqn 3.43 for effective r, but with the
                    % average correlation of the between-tree correlatons for those cores as r_bt and the
                    % actual number of cores as the effective number of cores per trees

                elseif nuse>1 & nuse==ntrees(n); % if more than one series in this year and the number of series equals number of 
                    % cores;  this means hae one core per tree, and should use r_bt as the effective
                    % correlation
                    % What rows, cols of the correlation matrix are needed?
                    J1=subfun01(find(L1));
                    irow=J1(:,1);
                    icol=J1(:,2);
                    ilinear=sub2ind(size(RR),irow,icol); % gives single-col index into correl matrix
                    r(n)=nanmean(RR(ilinear));
                    n2(n)= sum(~isnan(RR(ilinear))); % how many correlations go into the mean correlation
                    n3(n)=  nanmedian( NRR(ilinear)); % median sample length of pairs on which mean correl based
                elseif nuse>1; %  multiple series, more series than trees.  Use eqn 3.43 with the computed effective number of cores
                    % per tree (Ceff(n))
                    
                else;

                    r(n)=NaN;
                    n3(n)=NaN;
                    n2(n)=0;

                end;
                eflag=[];
            end; %  for n = 1:mX; % loop over years
            disp('here')


        elseif kmethod==3; % time invariant
            [mR,nR]=size(R);
            G=ones(mR,nR);
            L=logical(triu(G)) & Tdiff;
            Rsub=R(L); % cv of relevant correlations
            n2 = sum(~isnan(Rsub)); % number of correlations averaged for mean between-series r
            r = nanmean(Rsub);
            NRsub = NR(L); % cv of sample sizes for correlations
            n3 = nanmedian(NRsub); % median sample size for correlation
            eflag=[];

        end;
end; %  switch kwindow


%---------------------SUBFUNCTIONS

function J1=subfun01(j1)
%
% Build row and col indices to a reference matrix
% j1 are the col indices (variables) whose elements of matrix (row, col)
% are desired
% J1 are the rows and cols to pull from the correlation matrix

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
