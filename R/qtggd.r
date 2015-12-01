qtggd = function(p, logh=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE, res.approx=1e-4){
    mmax = logh + 2.5/b
    logm = seq(xmin, mmax, res.approx)
    cdf = ptggd(q=logm, logh=logh, a=a, b=b, xmin=xmin, lower.tail=lower.tail)
    icdf = approxfun(cdf, logm)
    p[p>1]=1
    p[p<0]=0
    if(log.p){p=exp(p)}
    return(icdf(p))
}