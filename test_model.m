% This model is based on Soto-Trevino et all 2005 and adds Q10's.
function test_model()
    AB = 1; % Enumerate cell types.
    PD = 2;
    Soma = 1; % Enumerate compartments.
    AIS  = 2; % (Axon initial segment.)
    
    % Model parameters %
    dt                  = 5*10^-3;  % [s] Time step
    sim_length          = 30.0;         % [s] Simulation length
    num_neurons         = 2;           % Number of neurons
    num_channels        = 10;          % Number of ion channel types.
    num_compartments    = 2;           % Number of compartments per neuron. Compartment 1 is soma. Compartment 2 is spike initiation zone.

    neurons(AB).compartments(AIS).C  =  1.5*10^-9; % [F] Capacitance of each compartment.
    neurons(AB).compartments(Soma).C =  9.0*10^-9; % [F]
    neurons(PD).compartments(AIS).C  =  6.0*10^-9; % [F]
    neurons(PD).compartments(Soma).C = 12.0*10^-9; % [F]
    
    % Leak conductances
    %              Soma  AIS
	g_leaks   = [  0.045 0.0018;   ...     % AB % [S] Leak conductances [AB-soma, AB-AIS; PD-soma PD-AIS] (Note: Compartment order different than for g_max_all.)
                   0.105 0.00081 ]*10^-6;  % PD
               
    % Leak reversal potentials
    %              Soma  AIS
    E_leak     = [ -50   -60;  ...         % AB [V]
                   -55   -55  ] * 10^-3;   % PD
    
    % Starting conditions
    starting_voltage = +20*10^-3; % [V] 
    for neuron = 1:num_neurons
        for compartment = 1:num_compartments 
            neurons(neuron).compartments(compartment).voltage = starting_voltage; 
        end
    end

    % History variables
    V_hist = zeros(num_neurons, num_compartments, ceil(sim_length/dt));
    I_hist = zeros(num_neurons, num_channels+1,   ceil(sim_length/dt)); % +1 is for leak channel.
    time_step = 1;
    
    %%%%%%%%%%%%%%%%%
    %   Main Loop   %
    %%%%%%%%%%%%%%%%%
    for sim_time = dt:dt:sim_length

        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find new voltage of each compartment %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        % Find currents within each neuron
        for neuron = 1:num_neurons
    
            ss_current      = zeros(1, num_compartments); % Sum of all currents at steady state  (active, leak, axial, gap and applied).
            conductance_sum = zeros(1, num_compartments); % Sum of all conductances (active, leak, axial and gap)
            
            for compartment = 1:num_compartments
                V = neurons(neuron).compartments(compartment).voltage;
                conductance_sum(compartment)  = g_leaks(neuron, compartment);
                leak_current =  g_leaks(neuron,  compartment)*(E_leak(neuron, compartment) -V);
                I_hist(neuron, 10+compartment, time_step) = leak_current;
                ss_current(compartment) = g_leaks(neuron,  compartment)*(E_leak(neuron, compartment)); 

                V_inf  =  ss_current(compartment)/conductance_sum(compartment);   % Weighted sum of reversal potentials.
                tau_V  =  neurons(neuron).compartments(compartment).C/(conductance_sum(compartment));   % Membrane time constant.
                neurons(neuron).compartments(compartment).voltage      =  V_inf + (V-V_inf)*exp(-dt/tau_V);
                %dV = leak_current/neurons(neuron).compartments(compartment).C*dt
                %neurons(neuron).compartments(compartment).voltage = V + dV;
                neuron
                V_inf
                ss_current(compartment)
                conductance_sum(compartment)
                V_hist(neuron, compartment, time_step) = neurons(neuron).compartments(compartment).voltage;
            end
        end
        
        time_step = time_step+1;
    end % End main loop    
    
    figure
    hold on
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(AB,Soma,:))*1000, 'r')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(AB,AIS, :))*1000, 'g')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(PD,Soma,:))*1000, 'b')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(PD,AIS, :))*1000, 'c')
        legend('AB Soma','AB AIS','PD Soma','PD AIS')
        xlabel('Time (ms)')
        ylabel('Membrane potential (mV)')
    hold off
    
    for neuron = 1:num_neurons
        figure
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron, 11,:))*1000, 'g.')
    end
end % End function