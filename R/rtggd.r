rtggd = function(n, logh=14, a=-1, b=1, xmin=10, res.approx=1e-4){
  if(length(n)>1){n=length(n)}
  return(qtggd(runif(n), logh=logh, a=a, b=b, xmin=xmin, res.approx=res.approx))
}
