function L= longrun(x)
% longrun: logical index to longest run of non-NaN in a vector
% function L = longrun(x)
% Last Revised 2006-3-23
%
% Logical index to longest run of non-NaN in a vector.  Needed by dtrendrw.m to pull
% longest stretch of OK core index from a series whose index blows up in places because
% of zero or negative growth curve. 
%
%****   INPUT 
%
% x(mX x 1)r vector of data, possibly with NaNs
%
%*** OUTPUT
%
% L (mx x 1)L vector of 1 or 0 telling whether x valid or NaN
%
%*** REFERENCES --- none
%
%*** TOOLBOXES NEEDED -- none
%
%*** UW FUNCTIONS CALLED -- none
%
%*** NOTES
%
% If a tie (2 stretches of same length NaN), earliest stretch is returned


if isscalar(x);
    error('x a scalar');
end
if size(x,2)~=1;
    error('x not col vector');
end;
if ~any(isnan(x));
    error('Do not call this function if no NaNs in vector');
end

%-- Attach a leading and trailing Nan to time series
xorig=x;
x = [NaN; x ; NaN];




% Find when series switches from NaN to OK and vice versa
x1 = x(1:(end-1));
x2 = x(2:end);
L1 = isnan(x1) & ~isnan(x2); % Mark years that are change from NaN to valid
L2 = ~isnan(x1) & isnan(x2); % mark years that are change from valid to NaN
i1 = find(L1); % index to elements of x marking chage from valid to NaN
i2 = find(L2); % index to elements of x marking chage from NaN to valid

% Find the longest stretch of non-NaN data
I1= repmat(i1',length(i2),1);
I2 = repmat(i2,1,length(i1));
d=(I2-I1);
d(d<=0)=NaN;
[s,s1]=nanmin(d);
[j,k]=max(s);

igo = i1(k); % start index of longest non-NaN stretch
isp = igo+j-1; % end index....

% Put results in logical vector
L=logical(zeros(length(xorig),1));
L(igo:isp)=1;




