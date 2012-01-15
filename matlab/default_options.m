function options = default_options()
% DEFAULT_OPTIONS
%   Return a struct with a default set of options.
%
options.n_moments = 2;        % two moments are fine
options.model_posterior = 1;  % and we want a model
			      % posterior
options.bprob = 1;            % we want break probabilities
options.utility = 1;          % compute utility gain
options.differential_entropy = 1;
options.multibin_entropy = 0;
options.effective_counts = 0;
options.marginal = 1;         % we also want the marginal

options.marginal_step = 0.01; % step size at which the
                              % density is evaluated
options.marginal_range(1) = 0.0;
options.marginal_range(2) = 1.0;

options.epsilon = 0.00001;    % precision for the extended
                              % prombs
options.threads = 1;          % no threads for the moment
options.stacksize = 256*1024; % some memory for prombs
options.algorithm = 0;        % 0: prombs, 1: no idea,
                              % 2: mgs
options.which = 0;            % compute everything for the
                              % first event (success)
options.samples(1) = 0;       % mgs burn in
options.samples(2) = 0;       % mgs samples

end % default_options
