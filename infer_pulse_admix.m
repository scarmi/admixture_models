function [best_t,best_m] = infer_pulse_admix(data,L,t_grid,m_grid)

% Infers the parameters of a pulse admixture model from the observed ancestry proportions.
% Usage:
% data: a vector with the observed *diploid* A ancestry proportion (i.e., the proorption of the two chromosomes that descend from the A
% population)
% L: the chromosome length (Morgans)
% t_grid: a list of values of t to search for. t is the time of admixture (generations ago)
% m_grid: a list of values of m to search for. m is the proportion of A ancestry contrinuted at the admixture event (the rest is B
% ancestry)

% The bin size for computing the likelihood of each data point (arbitrary)
dx = 0.01;
x = 0:dx:1;
h = histc(data,x);
h0 = sum(data==0);
% Subtract the number of zeros from the first bin
h(1) = h(1) - h0;
hall = [h0 h];

% Perform grid search
% Initialize log-likelihoods
lls = zeros(length(t_grid),length(m_grid));
for ti=1:length(t_grid)
    tt = t_grid(ti);
    for mi=1:length(m_grid)
        mm = m_grid(mi);
        [~,pdf] = pulse_diploid_ancestry_prop(dx,L,tt,mm);
        % Sum over all x grid points. For each point, count the log-likelihood for the number of times the observed data was within the bin
        ll = sum(hall .* log(pdf));
        lls(ti,mi) = ll;
    end
end
[best_t_ind,best_m_ind] = find(lls==max(lls(:)));
best_t = t_grid(best_t_ind);
best_m = m_grid(best_m_ind);
