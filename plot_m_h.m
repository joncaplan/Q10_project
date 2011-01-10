% Print out plots of steady state channel conductances
% (Note: code taken directly from model code and modified to plot steady state conditions.)

function plot_m_h
    min_voltage = -80; max_voltage = +30; voltage_step = 0.1; % [mV]

    ss_activations   = zeros(2,10,(max_voltage-min_voltage)/voltage_step  + 1);
    ss_inactivations = zeros(2,10,(max_voltage-min_voltage)/voltage_step  + 1);

    for neuron = 1:2
        for channel = 1:10

            index = 0;
            for V = min_voltage:0.1:max_voltage
                index = index+1;

                if neuron == 1, neuron_type = 1; else neuron_type = 2; end

                %Ca = Ca * 10^6; % Convert intracellular [Ca++] from Mol to μM. (Elsewhere [Ca++] is always in Mol.)
                Ca = 10000; % TEMP. FIX!!!!!!!!!!!!!!!!!!!!!!!!

                % Channel dynamics functions table.
                Na=1; K_AIS=2; CaT=3; CaS=4; nap=5; H=6; K_soma=7; KCa=8; A=9; proc=10; % Enumeration of all the currents.

                switch channel
                    case Na % I_Na (Axon initial segment)
                        m_inf = sigmoid(-1, 24.7, 5.29, V);
                        h_inf = sigmoid(+1, 48.9, 5.18, V);
                        tau_m = tau_sigmoid(1.32, -1.26, -1,  120, 25, V);
                        tau_h = tau_sigmoid(0,     0.67, -1, 62.9, 10, V)  * tau_sigmoid(1.5, 1, +1, 34.9, 3.6, V);

                    case CaT % I_CaT
                        m_inf = sigmoid(-1, 25, 7.2, V);
                        h_inf = sigmoid(+1, 36, 7  , V);
                        tau_m = tau_sigmoid(55, -49.5, -1, 58, 17, V);
                        switch neuron_type
                            case 1 % AB
                                tau_h = tau_sigmoid(87.5,  -75, -1, 50, 16.9, V);
                            case 2 % PD
                                tau_h = tau_sigmoid(350,  -300, -1, 50, 16.9, V);
                            otherwise 
                                STOP
                        end

                    case CaS % I_CaS
                        m_inf = sigmoid(-1, 22.5, 8.5, V);
                        tau_m = tau_sigmoid(16, -13.1, -1, 25.1, 26.4, V);

                    case nap % I_NaP
                        m_inf = sigmoid(-1, 26.8, 8.2, V);
                        h_inf = sigmoid( 1, 48.5, 4.8, V);
                        tau_m = tau_sigmoid(19.8, -10.7, -1, 26.5,  8.6,  V);
                        tau_h = tau_sigmoid(666,  -379,  -1, 33.6,  11.7, V);

                    case H % I_h
                        m_inf = sigmoid( 1, 70, 6, V);
                        tau_m = tau_sigmoid(272, +1499, -1, 42.2, 8.73, V);

                    case {K_AIS, K_soma} % I_K {axon initial segment, soma}
                        m_inf = sigmoid(-1, 14.2, 11.8, V);
                        tau_m = tau_sigmoid( 7.2, -6.4, -1, 28.3, 19.2, V);

                    case KCa % I_KCa
                        switch neuron_type
                            case 1
                                m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 4, V);
                            case 2
                                m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 8, V);
                            otherwise 
                                STOP
                        end        
                        tau_m = tau_sigmoid(90.3, -75.09, -1, 46, 22.7, V);

                    case A % I_A 
                        m_inf = sigmoid( -1, 27,   8.7,  V);
                        h_inf = sigmoid( +1, 56.9, 4.9,  V);
                        tau_m = tau_sigmoid( 11.6, -10.4, -1, 32.9, 15.2,  V);
                        tau_h = tau_sigmoid( 38.6, -29.2, -1, 38.9, 26.5,  V);    

                    case proc % I_proc
                        m_inf = sigmoid( -1, 12, 3.05,   V);
                        tau_m = 0.5*10^-3;

                    otherwise
                        STOP

                end
            ss_activations(neuron, channel, index)   = m_inf;
            ss_inactivations(neuron, channel, index) = h_inf;
            end
        end
    end

    for neuron = 1:2
        figure
        hold on

        plot(squeeze(ss_activations(neuron, :, :))')
        %legend('Na',  'K_A_I_S', 'CaT', 'CaS', 'nap', 'H', 'K_s_o_m_a')
        legend('Na',  'K_A_I_S', 'CaT', 'CaS', 'nap', 'H', 'K_s_o_m_a', 'KCa', 'A', 'proc')

        set(gca,'XTick', 1:10/voltage_step:(max_voltage-min_voltage)/voltage_step +1)
        voltages = ['-80';'-70';'-60';'-50';'-40';'-30';'-20';'-10';'  0';' 10';' 20'];        
        %voltages = ['-80';'-60';'-40';'-20';'  0';' 20'];
        set(gca,'XTickLabel',voltages)
        if neuron == 1 
            title('AB');
        else
            title('PD')
        end
    end
end


% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V)
    answer = 1/(1+exp(a*(V+b)/c));
end

function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = (a + b/(1+exp(c*(V+d)/e)))*10^-3;
end