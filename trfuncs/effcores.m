function [ntrees,ncores,ceff]=effcores(X,yrX,nms)
% effcores: effective number of cores each year
% [ntrees,ncores,ceff]=effcores(X,yrX,nms);
% Last revised 2006-3-20
%
% Effective number of cores each year.  Tree-ring time series matrix analyzed
% for number of trees, cores, and effective number of cores each year.  
%
%
%*** INPUT
%
% X(mX x nX)r time series
% yrX (mX x 1)i  time vector (e.g., year) for X
% nms (nX x 1)i integer "tree numbers" associated with each column of X (from treeid.m)
%
%*** OUTPUT
%
% ntrees (mX x 1)i  number of trees each year
% ncores (mX x 1)i number of cores each year
% ceff (mX x 1)r  effective number of cores each year (see Notes)
%
%*** REFERENCES -- none
%
% Briffa and Jones (1990): effective number of cores from eq 3.42, p.142
%
%*** UW FUNCTION CALLED
%
% trimnan -- trims leading and trailing all-NaN rows from X
% trailnan
%
%*** TOOLBOXES NEEDED
%
% Stats
%
%*** NOTES
%
% see References
%
% 2006-3-20 modified 


%--- CHECK INPUT


[mX,nX]=size(X);
ntemp=length(yrX);
if mX~=ntemp;
    error('yrX and X not same row size');
end;
yrr=yrX;


%---- GET THE BETWEEN-TREE INFORMATION

if isempty(nms);
    treexref=[];
else;
    treexref=nms;
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


%----- COMPUTE NUMBER OF SERIES IN EACH YEAR

ncores =   (sum(~isnan(X')))';




%-----  COMPUTE NUMBER OF TREES IN EACH YEAR


ntrees=zeros(mX,1); % initialize as no trees in any year

% Build matrix T1, same size as X, with tree number replacing data value
T1 = repmat(treexref',mX,1);
T1(isnan(X))=NaN;

% Find out how many trees in the file, and their arbitrary tree numbers
tree1 = unique(treexref); % col vector of arbitrary tree numers
ntree1 = length(tree1); % number of different trees in the data set

C=repmat(NaN,mX,ntree1); % to hold number of cores per year for each tree

% Loop over tree numbers, checking whether tree present in a year, and increment tree counter if so
for n =1:ntree1; % loop over unique arbitrary tree numbers
    thistree=tree1(n); % a tree number
    T2 = repmat(thistree,mX,nX);
    L = T2 == T1;
    L1=(any(L'))';
    ntrees = ntrees+L1;
    C(:,n)=(nansum(L'))';
end; % for n =1:ntree1; % loop over unique arbitrary tree numbers
clear T1 n thistree T2 L L1

% Compute effective number of cores per tree each year
C(C==0)=NaN; 
C1 = 1 ./ C; % 1 / ci, in notation of eq 3.42
C2=    (nansum(C1'))';
ntrees(ntrees==0)=NaN;
C2(C2==0) = NaN;
C3 = (1 ./ ntrees)  .* C2;
ceff = 1 ./ C3;



% figure(1);1
% hp1=plot(yrX,ncores,yrX,ntrees,yrX,Ceff);
% set(hp1(1),'Color',[.6 .6 1],'LineWidth',6);
% set(hp1(2),'Color',[1 0 0],'LineWidth',2);
% set(hp1(3),'Color',[0 0 0],'LineWidth',1);
% set(gca,'YLim',[0 1+max(ncores)]);
% ylabel('Number of series');
% xlabel('Time');
% legend('# cores','# trees','Ceff');
% grid on;
% title('Sample size in each year');
% [cL,cB,cW,cH]=figsize(.7,.5);
% set(gcf,'Position',[cL cB cW cH]);
