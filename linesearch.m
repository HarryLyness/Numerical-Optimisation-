function alpha = linesearch(x,fx,gx,s,theta,f)
% Backtracking line search algorithm
%
% input:    x - current iterate
%           fx - f(x)
%           gx - gradient f(x)
%           s - search direction (steepest descent is -gx)
%           theta - user threshold from Armijo condition
%           f - function handle
%
% output:   alpha

% Initialise alpha and set function increment for Armijo sufficient descent 
% condition
alpha = 1;
df = theta*gx'*s;

% Divide alpha by 2 until the sufficient descent condition is satisfied; 
% stop when alpha = 2^{-30} < 10^{-9}

while f(x+alpha*s) > fx + alpha*df 
%for i=1:30

    %if f(x+alpha*s) <= fx + alpha*df 
     %   break;
    %end
    
    alpha = alpha/2;
        
end

end

