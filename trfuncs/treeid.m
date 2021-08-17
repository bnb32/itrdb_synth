function [treenms,treexref,check]=treeid(nms)
% treeid: assign a unique tree "number" to each unique tree in rwl file
% function [treenms,treexref]=treeid(nms)
% Last Revised 2006-3-21
%
% Assign a unique tree "number" to each unique tree in rwl file.  Associate
% each input core id (nms) with a tree number.
% Called by member1.m as part of flowrec sequences.
%
%*** IN
%
% nms {n1 x 1}s  core ids in .rwl file
%
%*** OUT
%
% treenms{? x 1}s alphabetically ordered list of unique tree names, with site code attached
% treexref (n1 x 1)i cross-reference of each core in nms to a tree (row number) in treenms (Notes)
% check (? x ?)s table to check performance of program:
%
%*** REFERENCES --- none
%
%*** TOOLBOXES NEEDED -- none
%
%*** UW FUNCTIONS CALLED 
%
%*** NOTES
%
% User must follow a specific core-naming and core-numbering convention
% for treenum.m to work correctly. 
%
% treexref: the tree "number" in treexref is not necessarily the tree number in the core id.  For example,
% the file ut516.rwl, by Dave Grow, would be summarized as follows.  The core ids (nms) are listed at right.
% The assigned sequential tree number is at left.  Grow's tree #3 is sequential tree 2, and subsequent ids
% have different assigned tree numbers.  The higest number in treexref equals the number of trees in the
% dataset.  Thus even though Dave Grow lists trees up to "09", there are only 7 different trees represented.
%
% 1  CB01A
% 1  CB01B
% 2  CB03A
% 2  CB03B
% 3  CB04A
% 3  CB04B
% 4  CB05A
% 4  CB05B
% 5  CB06A
% 5  CB06B
% 6  CB07A
% 6  CB07B
% 7  CB09A
% 7  CB09B
% 
%
%
% CAUTION.  Problem if id has no core identifier.  For example,
%  pad12  would be converted to pad1, which assumes the 2 is a
%  core id.  In fact, pad12 might be intended as pad tree 12.
%  I haven't yet run across this as a problem.
%
% 2006-3-21 Modified to deal with Graybill AZ510.rwl, which has core names ending in two letters.
% Example
% SFP011AA 
% SFP012AA
% SPF021AA, etc
%  The trailing two letter appear to padding having nothing to  do with membership f the cores in trees.
%  Handle by substituting two blanks for all occurences of the trailing two letters
% Work with nms as cell



if isa(nms,'char');
    nms=cellstr(nms);
end;

ncores=length(nms);
treexref = repmat(NaN,ncores,1);
tnms = cell(ncores,1); % to hold tree name for each core in nms

%--- Deal with Graybill odd core-naming in which the core names all end with two letters. 
ctemp=char(nms); % cell to string matrix
L=isletter(ctemp); % mark chars in ctemp that are letters
Ltrail = L(:,(end-1):end); % are last two cols letters
if all(all(Ltrail)); % This is a Graybill two-letter
    ctemp(:,(end-1):end) = repmat('  ',ncores,1); % strip off trailing letters
    nms=cellstr(ctemp);
end;



for n = 1:ncores; % loop over cores
    
    nmthis =nms{n};
    charlast = nmthis(end);
	nmthis(end)=[]; % strip off what should be the core number or letter


	% Special case that had just stipped of a core-segment letter (e.g., from pad01a1)
    if ~isletter(charlast) & isletter(nmthis(end));
        nmthis(end)=[];
    end;
    if(length(nmthis)<3);
        tmp='XXX';
        tmp(1:length(nmthis))=nmthis;
        nmthis=tmp;
    end
	% Is there a letter in the remaining part of name other than chars 1-3
    a =nmthis;    
    
    
    a(1:3)=[];
    if any(isletter(a));
        warning(['Core ' nms{n} ' has letter after chars 1-3 after stripping off core letter']);
    end;
    tnms{n}=nmthis;
     
end; % for n = 1:ncores; % loop over cores

treenms=unique(tnms); % unique treenames
[tf,treexref] = ismember(tnms,treenms);
if ~all(tf);
    error('Some tree name not a member of the derived unique treenms; illogical');
end
check= [num2str(treexref) repmat('  ',length(nms),1)  char(nms)];
    