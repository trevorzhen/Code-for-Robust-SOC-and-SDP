% Author: Jianzhe Zhen
% Last Edited: 24-07-2020

function [A,b,nzer,nneg,npos] = core_fme(A,b,nvar)
% function [result] = fourmotz(A,b);
% Fourier-Motzkin Elimination for the problem A*x <= b
% with elimination of redundant inequalities
%   Input:  A ... matrix of size m x n
%           b ... vector of size m x 1
%           nvar (optional) ... number of variables that shall remain
%   Output: [A,b] reduced system of inequalities with all but nvar eliminated variables


if nargin < 3
	nvar = size(A,2)-1; 			% if nvar was not specified: choose 1
end

[m,n] = size(A);

if n > nvar
	temp = [A,b];			        % compound matrix [A|b]
    
	neg = find(temp(:,1)<0);	% find negative entries in first column
    nneg= size(neg,1);          % count negative entries in first column
    pos = find(temp(:,1)>0);	% find positive entries in first column
    npos= size(pos,1);          % count positive entries in first column
	zer = find(temp(:,1)==0);	% find zero entries in first column    
	nzer = sum(temp(:,1)==0);	% count zero entries in first column

	% divide according to entries in first column (see /Sch86a/: page 155, (16)):
	tempscalar=abs(temp([neg;pos],1)*ones(1,n+1));
    temp([neg;pos],:) =temp([neg;pos],:)./tempscalar;
	
    A=[];
	for i=1:npos
		A = [A;  ones(nneg,1)*temp(pos(i),2:n+1)+temp(neg,2:n+1)];
	end
	% append inequalities whose first column had a factor of zero:
	A = [A;temp(zer,2:n+1)];
	b = A(:,n);			     % extract vector b
	A = A(:,1:n-1);			 % reduce to the relevant A
	[A,b] = core_fme(A,b,nvar);		% recursive call
end