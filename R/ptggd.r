ptggd = function(q, logh=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE){
    p = gamma_inc((a+1)/b,10^(b*(q-logh)))/gamma_inc((a+1)/b,10^(b*(xmin-logh)))
    if(lower.tail){p=1-p}
    p[p>1]=1
    p[p<0]=0
    if(log.p){p=log(p)}
    return(p)
}