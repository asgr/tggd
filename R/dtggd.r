dtggd = function(x, logh=14, a=-1, b=1, xmin=10, log=FALSE){
  d = b*10^(a*(x-logh))*exp(-10^((x-logh)*b))/(10^logh * gamma_inc((a+1)/b,10^(b*(xmin-logh))))
  if(log){d=log(d)}
  return(d)
}
