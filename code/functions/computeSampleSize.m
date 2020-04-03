function [nout, tval, pval, power] = computeSampleSize(null_beta, beta, sigma2, gammaX, ncovariates, power, alpha)

% Two-sided t-test
tail = 0;
powerfun = @powerfunT;
significancefun = @significancefunT;

% Calculate one-sided Z value directly
sigma = sqrt(sigma2);
% out = z1testN(null_beta,beta,sigma,power,alpha,tail);

% Iterate upward from there for the other cases
N=ncovariates; %N = max(out,ncovariates); % t-test requires at least ncovariates
nout = searchupNextended(N,powerfun,significancefun,null_beta,beta,sigma,gammaX,ncovariates,power,alpha,tail);

df= nout-ncovariates;
tval = beta/sqrt(sigma2/df*gammaX);
pval = 2 * tcdf(-abs(tval), df); 

critL = tinv(alpha/2,df);   % note tinv() is negative
critU = -critL;
power = nctcdf(critL,df,tval) + nctcdf(-critU,df,-tval);


end


function power=powerfunT(mu0,mu1,sig,alpha,tail,n,ncovariates,gammaX)
%POWERFUNT T power calculation
    
    S = sig .* sqrt(gammaX./(n-ncovariates));       % std dev of mean
    ncp = (mu1-mu0) ./ S;     % noncentrality parameter
    df=n-ncovariates;
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

function pval=significancefunT(nout,beta,sigma,gammaX,ncovariates)
%SIGNIFICANCEFUNT T power calculation
    
df= nout-ncovariates;
tval = beta/(sigma*sqrt(gammaX/df));
pval = 2 * tcdf(-abs(tval), df); 
    
end

function N=searchupNextended(N,functP,functS,null_beta,beta,sigma,gammaX,ncovariates,desiredPower,alpha,tail)
%searchup Sample size calculation searching upward

    % Count upward until we get the value we need
    step_size = 2^7;
    todo = 0;
    while(~todo)
        N=N+step_size;
        actualpower = functP(null_beta,beta,sigma,alpha,tail,N,ncovariates,gammaX);
        actualSignificance = functS(N,beta,sigma,gammaX,ncovariates);
        todo = (actualpower > desiredPower) && (actualSignificance < alpha);
    end
    N=N-step_size;
    step_size=step_size/2;  
    
    for i_todo=1:7
        N=N+step_size;
        actualpower = functP(null_beta,beta,sigma,alpha,tail,N,ncovariates,gammaX);
        actualSignificance = functS(N,beta,sigma,gammaX,ncovariates);
        todo = (actualpower > desiredPower) && (actualSignificance < alpha);
        if todo
            N=N-step_size;
        end
        step_size=step_size/2;   
    end
    N=N+1;
end


%% 

function N=searchupN(N,F,mu0,mu1,sig,gammaX,ncovariates,desiredpower,alpha,tail)
%searchup Sample size calculation searching upward

    % Count upward until we get the value we need
    todo = 1:numel(alpha);
    while(~isempty(todo))
        actualpower = F(mu0,mu1(todo),sig,alpha(todo),tail,N(todo),ncovariates,gammaX);
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



function N=t1testN(mu0,mu1,sig,desiredpower,alpha,tail)
    %t1TESTN Sample size calculation for the one-sided Z test

    % Compute the one-sided normal value directly.  Note that we cannot do this
    % for the t distribution, because tinv depends on the unknown degrees of
    % freedom (n-1).
    if tail==0
        alpha = alpha/2;
    end
    todo=1;
    while(todo)
        actualSignificance = 0;
        todo = actualSignificance > desiredSignificance;
        N = N+1;
    end
   
    z1 = -norminv(alpha);
    z2 = norminv(1-desiredpower);
    mudiff = abs(mu0 - mu1) / sig;
    N = ceil(((z1-z2) ./ mudiff).^2);
end

