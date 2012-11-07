function result = lnchoose(n,k)
%
  result = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
%
end % lnchoose
