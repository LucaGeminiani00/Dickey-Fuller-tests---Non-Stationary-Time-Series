
%%%%%%%%%%%%%%%%%%%%%% Exercise 1 %%%%%%%%%%%%%%%%%%%%%%

AR1 = [1 -0.4];
rt = roots(fliplr(AR1))
poly(1./rt)
rt2 = roots(AR1)
poly(rt2)

%First method
ar = zeros(500,1); 
mu = 0;         % Mean
sigma2 = 0.2;     % Variance
epsilon = mu + sqrt(sigma2) * randn(500,1);
ar(1)=epsilon(1);

for t = 2:500
    ar(t) = 0.4*ar(t-1) +epsilon(t);
end 
figure; autocorr(ar,[],1)

%Second method
y1 = filter(1,AR1,epsilon)
figure; autocorr(y1,[],1)

%Comparison
y1 == ar;
if y1 == ar 
    disp('Identical Vectors')
else 
    disp('Some elements differ')
end 

%%%%%%%%%%%%%%%%%%%%%% Exercise 2 %%%%%%%%%%%%%%%%%%%%%%

%For loop
ar = zeros(500,1);
mu = 3;         % Mean
sigma2 = 0.4;     % Variance
epsilon= mu + sqrt(sigma2) * randn(500,1);
ar(1)=10;

for t = 2:500
    ar(t) = 0.6*ar(t-1)+epsilon(t);
end 

%Filter
AR1 = [1 -0.6];
epsilon(1)=10;
y1= filter(1,AR1,epsilon);

%%%%%%%%%%%%%%%%%%%%%% Exercise 3 %%%%%%%%%%%%%%%%%%%%%%

MA1 = [1 0.3];

%First method
ma = zeros(500,1); 
mu = 0;         % Mean
sigma2 = 0.3;     % Variance
epsilon = mu + sqrt(sigma2) * randn(500,1);
ma(1)=epsilon(1);

for t = 2:500
    ma(t) = 0.3*epsilon(t-1) + epsilon(t);
end 
figure; autocorr(ma,[],1)

%Second method
x1 = filter(MA1,1,epsilon);
figure; autocorr(x1,[],1)

%%%%%%%%%%%%%%%%%%%%%% Exercise 4 %%%%%%%%%%%%%%%%%%%%%%

%Function
%{
function[x,epsilon] = for_arma(T,sigma,phi,theta)

epsilon=sqrt(sigma) * randn(T, 1);
x=zeros(T,1);
x(1)=epsilon(1);
Q=length(theta);
P=length(phi);
A=zeros(P,1);
M=zeros(Q,1);

for t=2:T
    for p=1:P
        if t-1>=p
            A(p)=phi(p)*x(t-p);
        else
            A(p)=0;
        end
    end
    for q=1:Q
        if t-1>=q
            M(q)=theta(q)*epsilon(t-q);
        else
            M(q)=0;
        end
    end
    x(t)=sum(A)+sum(M)+epsilon(t);
end
end
%}
%Example 1
phi=[0.3 0.4];
theta=[0.3 0.4];
[x,epsilon] = for_arma(500,0.2,phi,theta);
figure;autocorr(x,[],2)

%Example 2
phi=[];
theta=[0.3 0.4 0.2];
[x,epsilon] = for_arma(500,0.2,phi,theta);
figure;autocorr(x,[],3)

%%%%%%%%%%%%%%%%%%%%%% Exercise 5 %%%%%%%%%%%%%%%%%%%%%%
%a)%Compute the empirical distributions of phi OLS using the data:  
%{
function [b,res,cov_b] = OLS(y,x)

b = inv(x'*x)*x'*y ; 
res = y-x*b;
cov_b = inv(x'*x)*cov(res);
end 
%}
%10000 simulations of AR(1) 250 OBS

emp_ols = zeros(10000,1); %vector of coefficients  
for i = 1:10000 
    ar = zeros(250,1);
    ar(1) = normrnd(0,0.4);
    epsilon = sqrt(0.4) * randn(250, 1);
    for t = 2:250
        ar(t) = 0.4*ar(t-1) + epsilon(t);
    end 
    arx = ar(1:249);  
    ary = ar(2:250); 
    [b,res,cov_b] = OLS(ary,arx);  
    emp_ols(i) = b ;  
end 

figure,histogram(emp_ols)

%b)(TEST WITH T = 200 ) 

emp_ols = zeros(10000,1); %vector of coefficients 
t_test = zeros(10000,1); %vector of t-test coefficients 
for g = 1:10000 
    ar = zeros(200,1);
    ar(1) = normrnd(0,0.2);
    epsilon = sqrt(0.2) * randn(200, 1);
    for t = 2:200
        ar(t) = 0.4*ar(t-1) + epsilon(t);
    end 
    arx = ar(1:199); %independent vars 
    ary = ar(2:200); %dependent vars
    [b,res,cov_b] = OLS(ary,arx); 
    emp_ols(g) = b;  
    t_test(g) = b/sqrt(cov_b) ;
end 

critical = tinv(1-(0.05/2),200-2); 
sum(abs(t_test)>critical)/10000 

%%%%%%%%%%%%%%%%%%%%%% Exercise 6 %%%%%%%%%%%%%%%%%%%%%%

G=[50,100,200,1000]; 
emp_ols_results = cell(length(G), 1);

for g= 1:length(G)
    T = G(g);
    emp_ols = zeros(10000,1); 
    for i = 1:10000 
        ar = zeros(T,1); 
        ar(1) = normrnd(0,0.2);
        epsilon = sqrt(0.2) * randn(T, 1);
        for t = 2:T
            ar(t) = 0.9*ar(t-1) + epsilon(t);
        end 
        arx = ar(1:T-1); %independent vars 
        ary = ar(2:T); %dependent vars
        [b,res] = OLS(ary,arx); 
        emp_ols(i) = b ;
    end 
    emp_ols_results{g} = emp_ols ; 
end 

%Comparison
for g = 1:length(G)
    T = G(g); 
    subplot(2, 2, g);
    histogram(emp_ols_results{g}, 30, 'Normalization','probability');
    title(['T = ' num2str(T)]);
end

%%%%%%%%%%%%%%%%%%%%%% Exercise 7 %%%%%%%%%%%%%%%%%%%%%%

%Empirical distribution OLS
a_ols = zeros(1000,1);
for i = 1:1000
    ma = zeros(250,1); 
    mu = 0;         % Mean
    sigma2 = 0.3;     % Variance
    epsilon = mu + sqrt(sigma2) * randn(250000,1);
    ma(1)=epsilon(1);
    for t = 2:250 
        ma(t) = epsilon(t) + 0.6*epsilon(t-1);
    end 
    xt_lagged = ma(1:249); %independent vars  
    xt = ma(2:250); %dependent vars
    [b,res] = OLS(xt,xt_lagged); 
    a_ols(i) = b ; 
end

%Mean
mean (a_ols)

%Increased observations
a_ols = zeros(1000,1);
for i = 1:1000
    ma = zeros(250000,1); 
    mu = 0;         % Mean
    sigma2 = 0.3;     % Variance
    epsilon = mu + sqrt(sigma2) * randn(250000,1);
    ma(1)=epsilon(1);
    for t = 2:250000 
        ma(t) = epsilon(t) + 0.6*epsilon(t-1);
    end 
    xt_lagged = ma(1:249999); %independent vars  
    xt = ma(2:250000); %dependent vars
    [b,res] = OLS(xt,xt_lagged); 
    a_ols(i) = b ; 
end
mean (a_ols)

%%%%%%%%%%%%%%%%%%%%%% Exercise 8 %%%%%%%%%%%%%%%%%%%%%%

%AR(1) process and the OLS sampling distribution  
emp_olsphi = zeros(1000,1); 
for i = 1:1000 
    ar = zeros(250,1);
    ar(1) = normrnd(0,0.2);
    epsilon = sqrt(0.2) * randn(250, 1);
    for t = 2:250
        ar(t) = ar(t-1) + epsilon(t);
    end 
    arx = ar(1:249);
    ary = ar(2:250) ; 
    [beta] = OLS(ary,arx); 
    emp_olsphi(i) = beta; 
end 

figure;histogram(emp_olsphi)

%t-test using a standard normal distribution
df_test = zeros(10000,1);

for i = 1:10000 
    ar = zeros(250,1);
    ar(1) = normrnd(0,0.2);
    epsilon = sqrt(0.2) * randn(250, 1);
    for t = 2:250
        ar(t) = ar(t-1) + epsilon(t);
    end 
    arx = [ones(249,1) ar(1:249)];
    ary = ar(2:250) - arx(:,2) ; %differenced dependent variable 
    [beta,res,cov_b] = OLS(ary,arx);
    df_test(i) = (beta(2)-0)/sqrt(cov_b(2,2));
end 

sum(df_test<-1.645)/10000
figure;histogram(df_test)

sum(df_test<-2.88)/10000

%comparison with DF distribution
percentiles = prctile(df_test, [1, 2.5, 5, 10]); 

%%%%%%%%%%%%%%%%%%%%%% Exercise 9 %%%%%%%%%%%%%%%%%%%%%%
%a) Build the Monte Carlo experiment; 

T= 250; %number of obs.
delta = 0.4; %drift parameter
for i = 1:10000 
    rw = zeros(T,1);
    rw(1) = normrnd(0,0.2);
    epsilon = sqrt(0.2) * randn(T, 1);
    for t = 2:T
        rw(t) = delta + rw(t-1) + epsilon(t);
    end 
    rwx = [ones(T-1,1) rw(1:T-1)];
    rwy = rw(2:T);
    [b,res,cov_b] = OLS(rwy,rwx); % b(1) is the est. of the intercept
    B_montecarlo(i,:) = b'; 
    df_test_1(i) = (b(2)-1)/sqrt(cov_b(2,2)); %the null is phi=1  
    df_test_2(i) = (b(1))/sqrt(cov_b(1,1));
end 

figure;histogram(B_montecarlo(:,1));
figure;histogram(B_montecarlo(:,2));
figure;histogram(df_test_1);
sum(df_test_1<-1.645)/10000

T= 1000; %number of obs.
delta = 0.4; %drift parameter
for i = 1:10000 
    rw = zeros(T,1);
    rw(1) = normrnd(0,0.2);
    epsilon = sqrt(0.2) * randn(T, 1);
    for t = 2:T
        rw(t) = delta + rw(t-1) + epsilon(t);
    end 
    rwx = [ones(T-1,1) rw(1:T-1)];
    rwy = rw(2:T);
    [b,res,cov_b] = OLS(rwy,rwx); % b(1) is the est. of the intercept
    B_montecarlo(i,:) = b'; 
    df_test_1(i) = (b(2)-1)/sqrt(cov_b(2,2)); %the null is phi=1  
    df_test_2(i) = (b(1))/sqrt(cov_b(1,1));
end 

figure;histogram(df_test_1);
sum(df_test_1<-1.645)/10000

%b) F-test
T = 250; 
phi = 1; 
alpha = 0.5; 
sigma = sqrt(0.2); 
ols_coeff_rho = zeros(1,1000); 
F_stat = zeros(1,1000); 
for i = 1:1000
    %AR(1) dataset 
    epsilon = sigma*randn(T,1); 
    y = zeros(T,1); 
    delta_y = zeros(T,1); 
    delta_y(1) = y(1); 
    for t = 2:T 
        y(t) = alpha + phi* y(t-1) + epsilon(t); 
        delta_y(t) = y(t) - y(t-1); 
    end 
    delta_y_t = delta_y(2:end); 
    lag = y(1:end-1); 
    mdl = fitlm(lag,delta_y_t); 
    SSR = mdl.SSE;  % Sum of squared residuals
    SST = mdl.SST;  % Total sum of squares
    df_R = mdl.NumCoefficients - 1;  % Degrees of freedom for the regression
    df_E = mdl.DFE;  % Degrees of freedom for the residuals
    F_stat(i) = ((SST - SSR) / df_R) / (SSR / df_E);
end 

critical = chi2inv(1-0.05, 1); 

sum(F_stat>critical)/1000 

df_chi2 = 1;

figure;
histogram(F_stat, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', 'blue', 'DisplayName', 'F-stat Histogram');
hold on;
x = 0:0.01:10; % Define a range of non-negative x values for the Chi-squared distribution
chi2_pdf = chi2pdf(x, df_chi2);
plot(x, chi2_pdf, 'r-', 'LineWidth', 2, 'DisplayName', 'Chi^2(1) Distribution');
%set x and y axis limits 
ylim([0, 1.5]);
xlim([0, max(x)]);
legend('Location', 'Best');
grid on;

hold off;

%c) deterministic time trend
simu = 1000;  
T = 250;  % Number of data points in each simulation

for j = 1:simu
    beta = 0.5;  % Set the time trend coefficient
    epsilon = randn(T, 1);  % Generate white noise
    
    % Create a time series with a deterministic time trend
    t = (1:T)';
    y = beta * t + cumsum(epsilon);
   
    % Perform the Dickey-Fuller test with the correct distribution
    [h, pValue,stat] = adftest(y, 'lags',1 , 'model', 'TS');
    emp_DF(j) = stat; 
end

histogram(emp_DF, 50, 'Normalization', 'probability')
xlabel('DF-stat')
ylabel('Probability')
title('Empirical Distribution DF-Stat')

percentiles = prctile(emp_DF, [1,2.5,5,10])

sum(emp_DF < -3.43)/simu