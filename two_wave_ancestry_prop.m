function [x,f] = two_wave_ancestry_prop(dx,L,t1,t2,m,mu)

% Comptues the PDF of x, the A ancestry proportion, for the two-wave admixture model
% Usage:
% dx: spacing between points (in the interval [0,1]) where the PDF is evaluated (see below for how x is generated). dx must satisfy
% 1/dx=integer.
% L: length of the chromosome (Morgans)
% t1: time of first admixture event (generations ago; A and B)
% t2: time of second admixture event (generations ago; A and the existing population)
% m: the proportion of A ancestry that enters the population at the time of the first admixture, t1 generations ago (the rest is B ancestry)
% mu: the proportion of A ancestry that enters the population at the time of the second admixture, t2 generations ago (the rest is the existing ancestry)
% Output:
% x: the grid of points where the PDF is evaluated, [0 (dx/2):(dx):(1-dx/2) 1]
% f: the PDF of the ancestry proportion
% Note that the PDF should satisfy sum(f)=1. The PDF at x=0 and x=1 is a delta function (mutiplied by the probability of no A ancestry or no
% B ancestry). For the remaining data points, the PDF is multiplied in dx, to generate the actual probability of an observation in the bin.


% Generate the vector of ancestry proportions where the PDF will be evaluated
x = [0 (dx/2):(dx):(1-dx/2) 1];
f = zeros(1,length(x));

% Generate a vector of La=xL.
Las = x*L;
% Evaluate the PDF at all internal data points
for i=2:(length(Las)-1)
    La = Las(i);
    f(i) = L*invlap('two_wave_laplace',L,[],1e-12,La,t1,t2,m,mu);
end

% Probability of A ancestry proportion = 0. Obtain directly from the probability to start at state B, and then stay in B for segment length
% longer than L.
pb_start = (1-m)*(1-mu);
pb_greater_than_L = exp(-L*(m*t1+mu*t2-m*mu*t2));
f(1) = pb_start * pb_greater_than_L;
% Same for the probability of A ancestry proportion = 1
pa_start = (mu + (1-mu)*m);
pa_greater_than_L = (-1).*exp(1).^((1/2).*L.*(((-1)+m).*t1+t2.*((-1)+mu+(-1).*m.*mu))).* ...
    ((-1)+m).*((-1)+mu).*((-1)+m+mu+(-1).*m.*mu).^(-1).*(t2.*mu+m.*(t1+( ...
    -1).*t2.*mu)).^(-1).*(4.*t1.*t2.*((-1)+m+mu+(-1).*m.*mu)+(t1+(-1).* ...
    m.*t1+t2+((-1)+m).*t2.*mu).^2).^(-1).*((m.*t1.*(((-1)+m).*t1+t2) ...
    .^2+(-1).*((-1)+m).*t2.*((1+3.*m.^2).*t1.^2+(-2).*t1.*t2+t2.^2).* ...
    mu+((-1)+m).^2.*((2+3.*m).*t1+(-2).*t2).*t2.^2.*mu.^2+(-1).*((-1)+m) ...
    .^3.*t2.^3.*mu.^3).*cosh((1/2).*L.*(4.*t1.*t2.*((-1)+m+mu+(-1).*m.* ...
    mu)+(t1+(-1).*m.*t1+t2+((-1)+m).*t2.*mu).^2).^(1/2))+(m.*t1.*(((-1)+ ...
    m).*t1+t2)+((-1)+m).*t2.*((-1).*(1+2.*m).*t1+t2).*mu+((-1)+m).^2.* ...
    t2.^2.*mu.^2).*(4.*t1.*t2.*((-1)+m+mu+(-1).*m.*mu)+(t1+(-1).*m.*t1+ ...
    t2+((-1)+m).*t2.*mu).^2).^(1/2).*sinh((1/2).*L.*(4.*t1.*t2.*((-1)+ ...
    m+mu+(-1).*m.*mu)+(t1+(-1).*m.*t1+t2+((-1)+m).*t2.*mu).^2).^(1/2)));
f(end) = pa_start * pa_greater_than_L;
f(2:(end-1)) = f(2:(end-1))*dx; % Easy to verify the distribution is normalized: sum(f)=1.
