function options = default_options()
% DEFAULT_OPTIONS
%   Return a struct with a default set of options.
%
options.n_moments = 2;        % two moments are fine
options.model_posterior = 1;  % and we want a model
			      % posterior
options.bprob = 1;            % we want break probabilities
options.kl_psi = 1;
options.kl_multibin = 0;
options.effective_counts = 0;
options.effective_posterior_counts = 0;
options.density = 1;          % we also want the density

options.density_step = 0.01;  % step size at which the
                              % density is evaluated
options.density_range(1) = 0.0;
options.density_range(2) = 1.0;

options.epsilon = 0.00001;    % precision for the extended
                              % prombs
options.threads = 1;          % no threads for the moment
options.stacksize = 256*1024; % some memory for prombs
options.algorithm = 0;        % 0: prombs, 1: mgs
options.which = 0;            % compute everything for the
                              % first event (success)
options.samples(1) = 0;       % mgs burn in
options.samples(2) = 0;       % mgs samples
options.hmm        = 0;       % do not use the hidden Markov model
options.rho        = 0.4;     % cohesion for the hidden Markov model

end % default_options
