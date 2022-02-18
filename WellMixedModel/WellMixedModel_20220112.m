
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate resistance evolution stochastically %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear
% close all

function output = WellMixedModel_20220112(adj_frac,pat)

%% Rates and Initialization

% Initial state 
% adj_frac = 0.1

% Rate constants (order: [sen res adj caf])
% Note: caf represents "carrying capacity" of CAFs, i.e. number of cells
% CAFs could have an effect on
ngr = [0.0238 0.0238 0.0238 0];   % net growth rate
d = [0.001 0.001 0.001 0];        % tonic death rate
b = ngr+d;                      % division rate
alpha_trt = [0.0410 0 0.0165 0];      % drug kill rate
alpha_pre = zeros(1,4);

mut = 7e-7;         % mutation rate [/division] (Loeb et al PNAS 2019)

rat_adj2sen = (adj_frac)/(1-adj_frac);
rec = ngr(1)*rat_adj2sen;

sen0=1e3;
pop0 = [sen0 0 sen0*adj_frac 1e5];
pop0 = [sen0 0 sen0*adj_frac sen0*adj_frac];
popmax = 1e9;
pop_treat = 1e9;

% rng(10)  % This controls the random number generator used in the simulation, so you get the same numbers every time you run the code (helpful for debugging)


% Initialize other variables (note that it's difficult to pre-allocate 
% variables to store the simulation outputs, because we don't know the
% size needed in a stochastic simulation)
pop = pop0;             % this variable stores the population at each iteration
pop_curr = pop(1,:);    % this variable keeps track of the current population structure

t = 0;                  % time [day]
t_curr = t;
tmax = 365*5;

% Initialize model variables
n = 1;          % iteration tracker
tau = NaN;      % initial tau
nSSAruns = NaN; % SSA run tracker (see description at line 62)
epsilon = .05;  % error tolerance in tau-leaping (used to calculate tau, see calculations in Cao Journal of Chemical Physics 2006 for more details)

treat = false;
alpha = alpha_pre;

%% Simulation    


while (t(end)<tmax) && (sum(pop(end,2:3))<popmax) && (sum(pop(end,:))>0) % run the simulation until 300 days have passed, the resistance population exceeds 1e3 cells, or the cells are eradicated
    
    % Identify the current state and time
    pop_curr = pop(end,:);
    t_curr = t(end);

    % Treatment conditions
    if sum(pop_curr(1:3))>=pop_treat && ~treat
        treat = true;
        t_trt = t_curr;
        n_trt = n;
        pop_trt = pop_curr;
        frc_adj = pop_trt(3)/(pop_trt(1)+pop_trt(3));
    end

    if treat
        alpha = alpha_trt;
    end
    
    
    % Define your propensity vector (describes the probability for each
    % event to occur; similar to your ODEs)
    evts = [b(1)*pop_curr(1); % sen division
            b(2)*pop_curr(2); % res division
            b(3)*pop_curr(3)*(1-pop_curr(3)/pop_curr(4)); % adj division --> adj
            b(3)*pop_curr(3)^2/pop_curr(4);   % adj division --> sen
            (d.*pop_curr)'; % natural death
            (alpha.*pop_curr)'; % drug kill
            mut*b(1)*pop_curr(1) % mutation
            rec*sum(pop_curr(1:3))]; % CAF recruitment
          
    a0 = sum(evts); % total sum of rates (describes how quickly events are occuring)
    
    % Define your state-change matrix (describes how the state/population
    % changes for each event; each row represents an event matching the
    % propensity vector evts, each column represents a population matching
    % the state vector pop, and each element in the matrix describes how
    % that population changes when that event occurs)
    numChange = [1 0 0 0; % division
                 0 1 0 0;
                 0 0 1 0;
                 1 0 0 0;  % spillover event
                 -eye(4); % natural death
                 -eye(4); % drug kill
                 0 1 0 0; % resistance mutation
                 0 0 0 1]; % CAF recruitment
             
    % Either run tau-leaping method or typical Gillespie SSA. Execute a
    % number of SSA runs if computed tau ever drops below expected time for
    % a single event (if this is the case, the tau-leaping method will run 
    % more slowly than Gillespie SSA). By default, it will start by trying
    % tau-leaping. See code at lines 97-104 that decide which method is 
    % appropriate.
    
    % Tau-leaping method (if Gillespie SSAs, code jumps to line 134)
    if isnan(nSSAruns)
        
        % Calculate tau, if not yet defined (by default, undefined)
        if isnan(tau)
            
            % Evaluate the Jacobian of your propensity vector (evts),
            % evaluated at the current state pop_curr. Note that in this
            % case, the Jacobian could have been defined outside the
            % simulation, since it doesn't change with respect to pop_curr.
            % That's not always the case though, so we'll leave it here.
            Jac = [diag(b(1:3)) zeros(3,1) % division
                   0 0 2*pop_curr(1)/pop_curr(4)*b(3) -pop_curr(1)^2/pop_curr(4)^2*b(3)
                   diag(d)
                   diag(alpha)
                   mut*b(1) 0 0 0
                   rec rec rec 0];
               
            % Calculate tau (read Cao Journal of Chemical Physics 2006
            % 'Efficient step size selection...' for description of tau
            % calculation.
            fjj = Jac*numChange';
            mu = fjj*evts;
            std = (fjj.^2)*evts;
            
            tau = min(min(epsilon*a0./(abs(mu)), epsilon^2*a0^2./(std)));
            
        end
        
        % Restart while loop iteration and switch to SSA if tau smaller 
        % than the expected time to next event (1/a0). Review Poisson
        % processes to understand why 1/a0 is the expected time to next
        % event.
        if tau < 1/a0
            nSSAruns = 10;  % counter for SSA runs
            tau = NaN;      % reset tau
            continue
            
        % Otherwise, calculate new state
        else
        
            % For each event, simulate number of times it occurs in time
            % interval tau and update the state
            pop_next = pop_curr;
            for i = 1:length(evts)
                num_evts = poissrnd(evts(i)*tau);
                pop_next = pop_next + numChange(i,:).*num_evts; % calculate change to state given that event 'i' occured 'num_evts' times
            end
            
            % Rerun with smaller tau if any populations are negative. This
            % is necessary, since tau-leaping is an 'approximate' and not
            % an 'exact' method. This means that the model could simulate
            % impossibilities, such as more death events than there are
            % cells in a population.
            if sum(pop_next<0)>0
                tau = tau/2;
                continue
            
            % Otherwise, save the new state
            else
                t_next = t_curr + tau;
                tau = NaN; % reset tau
            end
        
        end
        
    % Gillespie SSAs (no tau-leaping)
    else
        
        % Determine how long it takes for the next event to occur. This 
        % depends on the rate of events, and is drawn from an exponential 
        % distribution. Read up on Poisson processes to learn more. Read 
        % up on inverse transform sampling how delta_t is drawn.
        delta_t = -log(rand)/a0;
        t_next = t_curr + delta_t;
        
        % Determine which event in the propensity vector occurs (weighted by
        % their probabilities)
        evt_n = randsample(length(evts),1,true,evts);
        
        % How does that event change the state?
        pop_next = pop_curr + numChange(evt_n,:);
        
        % Reduce counter for SSA runs. If 0, reset counter (i.e. resume
        % tau-leaping method next iteration)
        nSSAruns = nSSAruns - 1;
        if nSSAruns < 1
            nSSAruns = NaN;
        end
        
    end
             
    % Update state and time variables
    n = n+1;
    pop(n,:) = pop_next;
    t(n) = t_next;
    % disp(log10(pop_curr))
    
end

output = [pat adj_frac t_trt pop_trt;
          pat adj_frac t_curr pop_curr];




%% Plot results
% 
% figure
% semilogy(t,pop(:,1:3))
% title('Resistance Evolution Simulation')
% xlabel('Time [days]')
% ylabel('Population Size [cells]')
% legend('Sensitive','Genetic Resistance','CAF-Mediated Resistance','Location','southeast')