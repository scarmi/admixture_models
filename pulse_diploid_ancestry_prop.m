function [x,y] = pulse_diploid_ancestry_prop(dx,L,t,m)

% Comptues the PDF of x, the *diploid* A ancestry proportion, for the pulse admixture model
% Usage:
% dx: spacing between points (in the interval [0,1]) where the PDF is evaluated (see below for how x is generated). dx must satisfy
% 1/dx=integer.
% L: the chromosome length (Morgans)
% t: the time of admixture (generations ago)
% m: the proportion of A ancestry contrinuted at the admixture event (the rest is B ancestry)
% Output:
% x: the grid of points where the PDF is evaluated, [0 (dx/2):(dx):(1-dx/2) 1]
% y: the PDF of the ancestry proportions
% Note that the PDF should satisfy sum(y)=1. The PDF at x=0 and x=1 is a delta function (mutiplied by the probability of no A ancestry or no
% B ancestry). For the remaining data points, the PDF is multiplied in dx, to generate the actual probability of an observation in the bin.


% Generate the vector of ancestry proportions where the PDF will be evaluated, except the boundaries (0 and 1)
x = (dx/2):(dx):(1-dx/2);

alpha = m*t;
lambda = (1-m)*t;

sq = sqrt(alpha*lambda*x.*(1-x));
% The haploid PDF of ancestry proportions
f = L*alpha*lambda/(alpha+lambda) * exp(-L*(lambda*x+alpha*(1-x))).*((alpha*x+lambda*(1-x)).*besseli(1,2*L*sq)./sq+2*besseli(0,2*L*sq));
f0 = lambda/(alpha+lambda)*exp(-alpha*L);
f1 = alpha/(alpha+lambda)*exp(-lambda*L);

n = length(x);
n2 = n/2; % Assuming n is even!!!
g = zeros(1,n);
% Implement convolution of two variables with PDF given by f, in order to get the PDF of their sum. The calculation is complicated by having 
% to compute separately the edge terms where x=0 or x=1
for i=1:n
    if i<=n2
        mn = 1;
        mx = 2*i-1;
        edge_term = f0 * (f(2*i-1)+f(2*i));
    else
        mn = 2*i-n;
        mx = n;
        edge_term = f1 * (f(2*i-n)+f(2*i-n-1));
    end   
    fs = f(mn:mx);
    g(i) = sum(fs .* fliplr(fs) * dx) + edge_term;
end
% The 2 is needed because the variable of interest is the average, not the
% sum, so when we change variables and divide by 2, we need to multiply g by 2.
g = g * 2*dx;
% This is the case when one individual is 0 and the other is 1, so the
% average is exactly 0.5.
g(n2+1) = g(n2+1) + 2*f0*f1;
% Probability that the diploid ancestry is 0 or 1.
g0 = f0^2;
g1 = f1^2;
% Output
x = [0 x 1];
y = [g0 g g1];