dtggd = function(x, scale=1e14, a=-1, b=1, xmin=1e10, log=FALSE){
  xtran=x/scale
  xmintran=xmin/scale
  if(log){
    d = log(b) + a*log(xtran) - xtran^b - log(scale) - log(gamma_inc((a+1)/b,xmintran^b))
  }else{
    d = b*(xtran^(a)*exp(-xtran^b))/(scale * gamma_inc((a+1)/b,xmintran^b))
  }
  return(d)
}

dtggd_log = function(x, scale=14, a=-1, b=1, xmin=10, log=FALSE){
  xtran=10^(x-scale)
  xmintran=10^(xmin-scale)
  if(log){
    d = log(log(10)) + log(b) + (a+1)*(x-scale)*log(10) - xtran^b - log(gamma_inc((a+1)/b,xmintran^b))
  }else{
    d = log(10)*b*(xtran^(a+1)*exp(-xtran^b))/gamma_inc((a+1)/b,xmintran^b)
  }
  return(d)
}

dtggd_ln = function(x, scale=log(1e14), a=-1, b=1, xmin=log(1e10), log=FALSE){
  xtran=exp(x-scale)
  xmintran=exp(xmin-scale)
  if(log){
    d = log(b) + (a+1)*(x-scale) - xtran^b - log(gamma_inc((a+1)/b,xmintran^b))
  }else{
    d = b*(xtran^(a+1)*exp(-xtran^b))/gamma_inc((a+1)/b,xmintran^b)
  }
  return(d)
}

ptggd = function(q, scale=1e14, a=-1, b=1, xmin=1e10, lower.tail=TRUE, log.p=FALSE){
    qtran=q/scale
    xmintran=xmin/scale
    p = gamma_inc((a+1)/b,qtran^b)/(gamma_inc((a+1)/b,xmintran^b))
    if(lower.tail){p=1-p}
    p[p>1]=1
    p[p<0]=0
    if(log.p){p=log(p)}
    return(p)
}

ptggd_log = function(q, scale=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE){
    qtran=10^(q-scale)
    xmintran=10^(xmin-scale)
    p = gamma_inc((a+1)/b,qtran^b)/gamma_inc((a+1)/b,xmintran^b)
    if(lower.tail){p=1-p}
    p[p>1]=1
    p[p<0]=0
    if(log.p){p=log(p)}
    return(p)
}

ptggd_ln = function(q, scale=log(1e14), a=-1, b=1, xmin=log(1e10), lower.tail=TRUE, log.p=FALSE){
    qtran=exp(q-scale)
    xmintran=exp(xmin-scale)
    p = gamma_inc((a+1)/b,qtran^b)/gamma_inc((a+1)/b,xmintran^b)
    if(lower.tail){p=1-p}
    p[p>1]=1
    p[p<0]=0
    if(log.p){p=log(p)}
    return(p)
}

qtggd = function(p, scale=1e14, a=-1, b=1, xmin=1e10, lower.tail=TRUE, log.p=FALSE, res.approx=1e-2){
    mmax = scale * 10^(2.5/b)
    logm = seq(log10(xmin), log10(mmax), res.approx)
    cdf = ptggd_log(q=logm, scale=log10(scale), a=a, b=b, xmin=log10(xmin), lower.tail=lower.tail)
    icdf = approxfun(cdf, logm)
    if(log.p){p=exp(p)}
    p[p>1]=1
    p[p<0]=0
    return(10^icdf(p))
}

qtggd_log = function(p, scale=14, a=-1, b=1, xmin=10, lower.tail=TRUE, log.p=FALSE, res.approx=1e-2){
    mmax = scale + 2.5/b
    logm = seq(xmin, mmax, res.approx)
    cdf = ptggd_log(q=logm, scale=scale, a=a, b=b, xmin=xmin, lower.tail=lower.tail)
    icdf = approxfun(cdf, logm)
    if(log.p){p=exp(p)}
    p[p>1]=1
    p[p<0]=0
    return(icdf(p))
}

qtggd_ln = function(p, scale=log(1e14), a=-1, b=1, xmin=log(1e10), lower.tail=TRUE, log.p=FALSE, res.approx=1e-2){
    mmax = scale + (2.5/b)*log(10)
    logm = seq(xmin/log(10), mmax/log(10), res.approx)
    cdf = ptggd_log(q=logm, scale=scale/log(10), a=a, b=b, xmin=xmin/log(10), lower.tail=lower.tail)
    icdf = approxfun(cdf, logm)
    if(log.p){p=exp(p)}
    p[p>1]=1
    p[p<0]=0
    return(icdf(p)*log(10))
}

rtggd = function(n, scale=1e14, a=-1, b=1, xmin=1e10, res.approx=1e-2){
  if(length(n)>1){n=length(n)}
  return(qtggd(runif(n), scale=scale, a=a, b=b, xmin=xmin, res.approx=res.approx))
}

rtggd_log = function(n, scale=14, a=-1, b=1, xmin=10, res.approx=1e-2){
  if(length(n)>1){n=length(n)}
  return(qtggd_log(runif(n), scale=scale, a=a, b=b, xmin=xmin, res.approx=res.approx))
}

rtggd_ln = function(n, scale=log(1e14), a=-1, b=1, xmin=log(1e10), res.approx=1e-2){
  if(length(n)>1){n=length(n)}
  return(qtggd_ln(runif(n), scale=scale, a=a, b=b, xmin=xmin, res.approx=res.approx))
}

tggd_mode=function(scale=1e14, a=-1, b=1,xmin=1e10){
  if(a<= 0){
    return(xmin)
  }else if(a==0){
    return(NaN)
  }else{
    return(scale*(a/b)^(1/b))
  }
}

tggd_mode_log=function(scale=14, a=-1, b=1,xmin=10){
  if(a<= -1){
    return(xmin)
  }else if(a==-1){
    return(NaN)
  }else{
    return(scale+log10((a+1)/b)^(1/b))
  }
}

tggd_mode_ln=function(scale=log(1e14), a=-1, b=1,xmin=log(1e10)){
  if(a<= -1){
    return(xmin)
  }else if(a==-1){
    return(NaN)
  }else{
    return(scale+log((a+1)/b)^(1/b))
  }
}

tggd_mean=function(scale=1e14, a=-1, b=1, xmin=1e10){
  xmintran=xmin/scale
  return(scale*gamma_inc((a+2)/b,xmintran^b)/gamma_inc((a+1)/b,xmintran^b))
}
