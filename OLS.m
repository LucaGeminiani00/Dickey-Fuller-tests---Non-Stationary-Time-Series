function [b,res,cov_b] = OLS(y,x)

b = inv(x'*x)*x'*y ; 
res = y-x*b;
cov_b = inv(x'*x)*cov(res);

end 