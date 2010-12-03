% This model is based on Soto-Trevino et all 2005 and adds Q10's.

function model()
    AB = 1; % Enumerate cell types.
    PD = 2;
    Soma = 1; % Enumerate compartments.
    AIS  = 2; % (Axon initial segment.)
    
    % Model parameters %
    dt                  = 0.05;    % [ms] Time step
    sim_length          = 10*1000; % [ms] Simulation length
    num_neurons         = 2;       % Number of neurons
    num_channels        = 10;      % Number of ion channel types.
    compartment_num     = 2;       % Number of compartments per neuron. Compartment 1 is soma. Compartment 2 is spike initiation zone.
    Ca_steady_state     = 0.5;     % [uM]
    g_gap               = 0.75;    % [uS] Conductance between the soma of the two neurons.
    neurons(AB).g_axial = 0.3;     % [uS] Conductance between the two compartments of the AB neuron.
    neurons(PD).g_axial = 1.05;    % [uS] Conductance between the two compartments of the PD neuron.
    
    neurons(AB).compartments(AIS).C  =  1.5; % [nF] Capacitance of each compartment.
    neurons(AB).compartments(Soma).C =  9.0; % [nF]
    neurons(PD).compartments(AIS).C  =  6.0; % [nF]
    neurons(PD).compartments(Soma).C = 12.0; % [nF]
    
    
    R = 8.31447215;  % [J/mol/K]
    T = 273;         % [Kelvin] %%%%%%%%%%%%% FIX
    z = +2;          % Charge of Ca++.
    F = 96485.3399;  % [C/mol]

    
    E_Na   =  50; % [mV] % Reversal potentials.
    E_K    = -80;
    E_Ca   =   0; % (Placeholder. Nernst calculation later.)
    E_Nap  =  50;
    E_H    = -20;
    E_KCa  = -80;
    E_A    = -80;
    E_proc =   0;

    E_all = [E_Na E_K E_Ca E_Ca E_Nap E_H E_K E_KCa E_A E_proc]; % Store reversal potentials in vector which matches the channel vector. (Ca values are placeholders.)

    g_max_all = [ 300  52.5      55.2  9   2.7   0.054  1890    6000    200    570; ... % AB % [uS] Conductance values for all (non-leak) channels.
                 1100  150       22.5  60  4.38  0.219  1576.8  251.85  39.42  0   ];   % PD
             
	g_leaks   = [  0.045 0.0018;   ... % AB % [uS] Leak conductances [AB-soma, AB-AIS; PD-soma PD-AIS] (Note: Compartment order different than for g_max_all.)
                   0.105 0.00081 ];    % PD
               
    
    
    % Set up model starting conditions
    starting_voltage = -70; % [mV] 
    for neuron = 1:num_neurons
        for compartment = 1:compartment_num 
            neurons(neuron).compartments(compartment).voltage = starting_voltage;
        end
    end
    
    for neuron = 1:num_neurons
        neurons(neuron).Ca = Ca_steady_state;
    end
    

    %%%%%%%%%%%%%%%%%
    %   Main Loop   %
    %%%%%%%%%%%%%%%%%
    for sim_time = dt:dt:sim_length

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find channel activation and inactivations  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for neuron = 1:num_neurons
            for channel = 1:num_channels
                if channel <= 2, compartment = 2; else compartment = 1; end % First two channels are in axon initial segment.
                [m h] = get_channel_state(neuron, channel, compartment);
                neurons(neuron).channels(channel).m = m;
                neurons(neuron).channels(channel).h = h;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % Calculate currents between compartments and neurons. %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        % First get the voltage differences between compartments and neurons.
        for neuron = 1:num_neurons
            delta_V_compartments(neuron) = neurons(neuon).compartments(1).voltage - neurons(neuon).compartments(2).voltage;
        end
        delta_V_gap = neurons(1).compartments(1).voltage  - neurons(2).compartments(1).voltage; % Compare somatic (compartment 1, where gap juntions are located) voltages.
        
        % Then calculate the currents.
        % ... within neuron
        for neuron = 1:num_neurons
            neurons(neuron).I = delta_V_compartments(neuron) * neurons(neuron).g_axial;
        end
        % ... between neurons.
        I_gap = delta_V_gap * g_gap;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate channel currents  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for neuron = 1:num_neurons
            % Find Ca reversal potential via Nernst equation.
            E_Ca = R*T/(z*F)*log(Ca_out/Ca_in);
            E_all(3) = E_Ca; E_all(4) = E_Ca; % Set reversal values for Ca++.
            
            for channel = 1:num_channels
                m = neurons(neuron).channels(channel).m;
                h = neurons(neuron).channels(channel).h;
                if (channel<=2), V = neurons(neuron).compartments(2).voltage; else V = neurons(neuron).compartments(1).voltage; end % Get voltage from correct compartment.
                g_max = g_max_all(1,1); % FIX
                switch channel % Lookup the channel activation and inactivation exponents.
                    case 1      % I_Na
                        a = 3; b = 1; 
                    case 2      % I_K (AIS)
                        a = 4; b = 0;
                    case 3      % I_CaT
                        a = 3; b = 1;
                    case 4      % I_CaS
                        a = 3; b = 0;
                    case 5      % I_nap
                        a = 3; b = 1;
                    case 6      % I_h
                        a = 1; b = 0;
                    case 7      % I_k (soma)
                        a = 4; b = 0;
                    case 8      % I_KCa
                        a = 4; b = 0;
                    case 9      % I_A
                        if (neuron == AB), a = 3; elseif (neuron == PD), a = 4; else STOP; end
                               b = 1;
                    case 10     % I_proc
                        a = 1; b = 0;
                end
                        
                g_channel = g_max * m^a * h^b;
                neurons(neuron).channels(channel).I = g_channel * (V - E_all(channel));
            end
            for compartment = 1:2 % Calculate leak currents
                I_leak(neuron, compartment) = g_leaks(neuron, compartment) * neurons(neuron).compartments(compartment).voltage;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find new voltage of each compartment %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        for neuron = 1:num_neurons
            total_current = axial_current + coupling_current + channel_current;
            %dV = current and capacitance;
            neurons(neuron).compartments(1).voltage = neurons(neuron).compartments(1).voltage + dV;
            neurons(neuron).compartments(2).voltage = neurons(neuron).compartments(1).voltage + dV;
        end
        
    end
end


%%%%%%%%%%%%%%
% Functions  %
%%%%%%%%%%%%%%

function [m, h] = get_channel_state(neuron, channel, compartment)
    V = neurons(neuron).voltage(compartment); % [mV]
    if neuron == 1, neuron_type = 1; else neuron_type = 2; end

    switch channel
        case 1 % I_Na (Axon initial segment)
            m_inf = sigmoid(-1, 24.7, 5.29, V);
            h_inf = sigmoid(+1, 48.9, 5.18, V);
            tau_m = tau_sigmoid(1.32, -1.26, -1,  120, 25, V);
            tau_h = tau_sigmoid(0,     0.67, -1, 62.9, 10, V)  * tau_sigmoid(1.5, 1, +1, 34.9, 3.6, V);

        case 3 % I_CaT
            m_inf = sigmoid(-1, 25, 7.2, V);
            h_inf = sigmoid(+1, 36, 7  , V);
            tau_m = tau_sigmoid(55, -49.5, -1, 58, 17, V);
            switch neuron_type
                case 1
                    tau_h = tau_sigmoid(87.5,  -75, -1, 50, 16.9, V);
                case 2
                    tau_h = tau_sigmoid(350,  -300, -1, 50, 16.9, V);
                otherwise 
                    STOP
            end
            
        case 4 % I_CaS
            m_inf = sigmoid(-1, 22.5, 8.5, V);
            h_inf = 1;
            tau_m = tau_sigmoid(16, -13.1, -1, 25.1, 26.4, V);
            tau_h = 1;
            
        case 5 % I_NaP
            m_inf = sigmoid(-1, 26.8, 8.2, V);
            h_inf = sigmoid( 1, 48.5, 4.8, V);
            tau_m = tau_sigmoid(19.8, -10.7, -1, 26.5,  8.6,  V);
            tau_h = tau_sigmoid(666,  -379,  -1, 33.6,  11.7, V);
            
        case 6 % I_h
            m_inf = sigmoid( 1, 70, 6, V);
            h_inf = 1;
            tau_m = tau_sigmoid(272, +1499, -1, 42.2, 8.73, V);
            tau_h = 1;
            
        case {2, 7} % I_K {axon initial segment, soma}
            m_inf = sigmoid(-1, 14.2, 11.8, V);
            h_inf = 1;
            tau_m = tau_sigmoid( 7.2, -6.4, -1, 28.3, 19.2, V);
            tau_h = 1;
            
        case 8 % I_KCa
            Ca = neurons(neuron).Ca;
            switch neuron_type
                case 1
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 4);
                case 2
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 8);
                otherwise 
                    STOP
            end        
            h_inf = 1;
            tau_m = tau_sigmoid(90.3, -75.09, -1, 46, 22.7, V);
            tau_h = 1;
            
        case 9 % I_A 
            m_inf = sigmoid( -1, 27,   8.7,  V);
            h_inf = sigmoid( +1, 56.9, 4.9,  V);
            tau_m = tau_sigmoid( 11.6, -10.4, -1, 32.9, 15.2,  V);
            tau_h = tau_sigmoid( 38.6, -29.2, -1, 38.9, 26.5,  V);    
            
        case 10 % I_proc
            m_inf = sigmoid( -1, 12, 3.05,   V);
            h_inf = 1;
            tau_m = 0.5;
            tau_h = 1;
            
        otherwise
            STOP
    end
    % Integrate to find m and h at next time step. 
    m = m_inf * tau_m; %%%%%%%%% Placeholder code %%%%%%%%%%
    h = h_inf * tau_h;
end

% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V)
    answer = 1/(1+exp(a*(V+b)/c));
end

function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = a + b(1+exp(c*(V+d)/e));
end