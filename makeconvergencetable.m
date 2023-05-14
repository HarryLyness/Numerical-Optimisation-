function [conv_table,e] = makeconvergencetable(xex,xsol,nsol)
% function to make the testing convergence tables
%
% Inputs: 
%         xex - exact solution 
%         xsol - computed iterates x_n
%         nsol - number of iterates x_0,...,x_n
%
% Outputs: 
%         conv_table - convergence table 
%         e - error vector
%
% calculates the error at each x_n
for i=1:nsol+1
    e(i) = norm(xsol(:,i)-xex);
end

% is the error quadratic? computes e_k+1/e_k^2
quadconvergence = zeros(1,nsol+1-1);
quadconvergence(1:nsol) = e(2:nsol+1)./(e(1:nsol).^2);
quadconvergence = quadconvergence';

% is the error superlinear? computes e_k+1/e_k
superlinear = zeros(1,nsol+1-1);
superlinear(1:nsol) = e(2:nsol+1)./(e(1:nsol));
superlinear = superlinear';

% outputs table, note that outputs n-1 iterations for table 
% since the n+1 iterate is not computed... 
n = [0:nsol-1]';
conv_table = table(n,quadconvergence,superlinear);
%re-lableing headers of table
conv_table = renamevars(conv_table,['n'],['iteration']);
conv_table = renamevars(conv_table,['quadconvergence'],['||e(k+1)||_2/||e(k)||_2^2']);
conv_table = renamevars(conv_table,['superlinear'],['||e(k+1)||_2/||e(k)||_2']);
end