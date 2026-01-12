function [spike_times,spike_ids,Irecons] = TauLeaping_Adaptive_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state,params)
% This is a basic implementation of the Tau Leaping Algorithm.
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
eps = 0.03;
nc = 10;

N = size(W,1);
Tau = 0.1;
time_points = t_min:Tau:t_max;
active = init_state(1,:);
network_states = zeros(N, length(time_points));
event_count = zeros(1, length(time_points));

propensity = beta .* active .* feval(response_fn, W*active' + I(:,1))'
 + alpha .* ~active;

network_states(:,1) = active;
event_count(1) = sum(active == 1);
for t in time_points
    Icurr = I(:,floor(t+1));
    
    % Calculate propensities
    propensity = beta .* active .* feval(response_fn, W*active' + Icurr)' ...
        + alpha .* ~active;
    
    % Adaptive Tau selection
    if any(num_events > 0)
        Lj = argmin(state_change / propensity);
        criticality = (Lj <= nc) && (propensity > 0);
        proder = sum(propensity/states * state_change)
        mu = 0
        var2 = 0 
        taup = min([0, 0])
        taupp = poissrnd(1/sum(propensity));

        if taup*sum(propensity) <= 1 
            print('Tau Leaping is not efficient with the current dynamics.\n Consider using Gillespie Algorithm instead.\n');
            return;
        end
    end

    if taup < taupp
        Tau = taup;
        num_events = poissrnd(propensity * Tau);
        num_events(criticality) = 0;
    elseif taup >= taupp
        Tau = taupp;
        num_events = poissrnd(propensity * Tau);
        j = find(propensity ./ sum(propensity) > rand, 1, 'first');
        num_events(j) = 1;
    end

    % Update states
    for i = 1:N
        if num_events(i) > 0 
            if active(i) == 0
                active(i) = num_events(i) // 2; % neuron becomes active
                propensity = propensity + num_events(i) * (W(:,i)' .* beta .* feval(response_fn, W*active' + Icurr)');
            else
                active(i) = num_events(i) // 2; % neuron becomes quiescent
                propensity = propensity - num_events(i) * (W(:,i)' .* beta .* feval(response_fn, W*active' + Icurr)');
            end
        else
            print('No event occurred for neuron %d between time %.2f and %.2f\n', i, t, t+Tau');
        end
        propensity = beta .* active .* feval(response_fn, W*active' + I(:,1))'
            + alpha .* ~active;
    end
    network_states(:,t) = active;
    event_count(t) = sum(active == 1);
end

spike_times = event_count .* time_points;
spike_ids = network_states;
end