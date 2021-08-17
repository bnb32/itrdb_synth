function z = temp(X, varargin)
[M,N] = size(X);

if nargin ==1
    nf =N*4;
    norm = 'Y';
elseif nargin ==2  
    nf = varargin{1};
    norm = 'Y';
elseif nargin ==3
    if isempty(varargin{1})
        nf = N*4;
    else
        nf = varargin{1};
    end
    norm = varargin{2};
else
    error('Invalid number of input arguments')
end
disp(nf);
disp(norm);
z = nf;