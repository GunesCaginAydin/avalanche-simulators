function [spike_times,spike_ids,Irecons] = TauLeaping_Adaptive_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state,params,mode)
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
active = init_state(1,:)';
network_states = zeros(N, length(time_points));
event_count = zeros(N, length(time_points));

network_states(:,1) = active;
event_count(:,1) = zeros(N,1);
count = 1;

current = W*active + I(:,1);
propensity = beta .* ~active .* feval(response_fn, current) ...
    + alpha .* active;
    
for t in time_points
    num_events = poissrnd(propensity * Tau);
    
    if any(num_events > 0)
        Lj = argmin(state_change / propensity);
        criticality = (Lj <= nc) && (propensity > 0);
        proder = sum(propensity/states * state_change)
        mu = 0
        var2 = 0 
        taup = min([0, 0])
        taupp = poissrnd(1/sum(propensity));

        if taup*sum(propensity) <= 1 
            print("Tau Leaping is not efficient with the current dynamics.\n
                        Consider using Gillespie's Algorithm instead.\n");
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

    for i = 1:N
        if num_events(i) > 0 
            if (num_events(i) // 2)
                if active(i) == 1
                    current += W(:,i);
                    active(i) = 0;
                else
                    current -= W(:,i);
                    active(i) = 1;
                end
            else
                continue;
            end
        else
            print('No event occurred for neuron %d between time %.2f and %.2f\n', i, t, t+Tau');
        end
    end
    count += 1;

    Ipast = I(:,floor(t+1));
    Icurr = I(:,floor(t+Tau+1));

    current = W*active + Icurr - Ipast;

    propensity = beta .* ~active .* feval(response_fn, current) ...
        + alpha .* active;
    network_states(:,count) = active;
    event_count(:,count) = num_events;
end

if mode == 0
    print('Accumulating spike times at the end of each time leap.\n')
    time_points = time_points
elseif mode == 1
    print('Uniformly distributing spike times throughout each time leap.\n')
    trecons = [];
    for i=1:size(event_count,2)
        trecons = [trecons linspace(time_points(i), time_points(i)+Tau, sum(event_count(:,i)))];
    end
    time_points = trecons
elseif mode == 2
    print('Uniformly distributing spike times with Gaussian Noise throughout each time leap.\n')
    trecons = [];
    for i=1:size(event_count,2)
        trecons = [trecons linspace(time_points(i), time_points(i)+Tau, sum(event_count(:,i))...
            + rand(1, sum(event_count(:,i))*0.1*Tau))]
        if trecons[1]<=0
            trecons[1] = 0;
        end
    end
    time_points = trecons
end

spike_times = [];
for i=1:len(time_points)
    spike_times = [spike_times repmat(time_points(i), 1, sum(event_count(:,i)))];
end
spike_ids = network_states;
end