function [nout, tval, pval, power] = computeSampleSize(null_value, standard_error, value, power, alpha, ncovariates)

% Two-sided t-test
tail = 0;
powerfun = @powerfunT;

% Calculate one-sided Z value directly
out = z1testN(null_value,value,standard_error,power,alpha,tail);

% Iterate upward from there for the other cases
out = max(out,2); % t-test requires at least 2
nout = searchupN(out,powerfun,null_value,value,standard_error,power,alpha,tail,ncovariates);

df= nout-ncovariates+1;
tval = value / (standard_error./ sqrt(nout));
pval = 2 * tcdf(-abs(tval), df); 

critL = tinv(alpha/2,df);   % note tinv() is negative
critU = -critL;
power = nctcdf(critL,df,tval) + nctcdf(-critU,df,-tval);


end


function power=powerfunT(mu0,mu1,sig,alpha,tail,n,ncovariates)
%POWERFUNT T power calculation
    
    S = sig ./ sqrt(n);       % std dev of mean
    ncp = (mu1-mu0) ./ S;     % noncentrality parameter
    df=n-ncovariates+1;
    if tail==0
        critL = tinv(alpha/2,df);   % note tinv() is negative
        critU = -critL;
        power = nctcdf(critL,df,ncp) + nctcdf(-critU,df,-ncp); % P(t < critL) + P(t > critU)
    
    elseif tail==1
        crit = tinv(1-alpha,df);
        power = nctcdf(-crit,df,-ncp); % 1-nctcdf(crit,n-1,ncp), P(t > crit)
    
    else % tail==-1
        crit = tinv(alpha,df);
        power = nctcdf(crit,df,ncp); % P(t < crit)
    end        
end

function N=searchupN(N,F,mu0,mu1,args,desiredpower,alpha,tail, ncovariates)
%searchup Sample size calculation searching upward

    % Count upward until we get the value we need
    todo = 1:numel(alpha);
    while(~isempty(todo))
        actualpower = F(mu0,mu1(todo),args,alpha(todo),tail,N(todo), ncovariates);
        todo = todo(actualpower < desiredpower(todo));
        N(todo) = N(todo)+1;
    end
end

function N=z1testN(mu0,mu1,sig,desiredpower,alpha,tail)
    %Z1TESTN Sample size calculation for the one-sided Z test

    % Compute the one-sided normal value directly.  Note that we cannot do this
    % for the t distribution, because tinv depends on the unknown degrees of
    % freedom (n-1).
    if tail==0
        alpha = alpha/2;
    end
    z1 = -norminv(alpha);
    z2 = norminv(1-desiredpower);
    mudiff = abs(mu0 - mu1) / sig;
    N = ceil(((z1-z2) ./ mudiff).^2);
end


