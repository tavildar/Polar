
%main_de


% close all;
% clear all;
addpath('GaussianApproximation/');
addpath('CapacityHelper/');

initialize_phi(10^(-5));
n = 10;
block_length = 2^n;
tic

target_bler = 10^(-5);
rate_vec = [1/32, 1/16, 1/8, 1/4, 2/4, 3/4, 7/8];
snr_increment = 0.25;
snr_db_vec = -10:snr_increment:30;

constellation_vec = {'ask4-gray', 'ask4-sp', 'ask16-gray', 'ask16-sp'}; 
receiver_vec = {'bicm', 'mlc','bicm', 'mlc'}; 

snr_needed = zeros(length(rate_vec), length(constellation_vec));
ebno_needed = zeros(length(rate_vec), length(constellation_vec));

figure(2);
for i_const = 1 : length(constellation_vec)

    constellation_name = constellation_vec{i_const};
    
    reciver_algo = receiver_vec{i_const};
    
    modulation = Constellation(constellation_name);
    
    starting_snr_index = 1;

    for rate_index = 1 : length(rate_vec)
        rate = rate_vec(rate_index);
        
        disp(['Processing constellation : ', constellation_name, ' at rate : ', num2str(rate)]);
        polar_code = PolarCode(block_length, ceil(rate*block_length), 0.5);
        
        for snr_index = starting_snr_index : length(snr_db_vec)
            
            design_snr_db = snr_db_vec(snr_index);
            bler_estimate = polar_code.ga_code_construction(design_snr_db, constellation_name, reciver_algo);
            if bler_estimate < target_bler
                break;
            else
                prev_bler_estimate = bler_estimate;
            end
           
        end
        
        if snr_index == starting_snr_index
            disp('Possibly too high starting SNR');
            k = waitforbuttonpress;
        end

        starting_snr_index = max(snr_index - 1, 1);
        
        snr_needed(rate_index, i_const) = (snr_db_vec(snr_index) * log(prev_bler_estimate/target_bler) + snr_db_vec(snr_index - 1) * log(target_bler/bler_estimate) ) / log(prev_bler_estimate/bler_estimate) ;
        ebno_needed(rate_index, i_const) = snr_needed(rate_index, i_const)  - 10*log10(rate) - 10*log10(modulation.n_bits);

    end
    
    plot(ebno_needed(:, i_const), rate_vec*modulation.n_bits, 'LineWidth', 2);
    hold on; grid on;
    legend_vec{i_const} = [constellation_vec{i_const}, ' - ', receiver_vec{i_const}];

end
legend(legend_vec);
xlabel('EbNo (dB)', 'FontSize', 14);
ylabel('Rate', 'FontSize', 14);
title(['GA E_b/N_0 vs Rate performance for target bler = ', num2str(target_bler)], 'FontSize', 14);
savefig('gauss-approx.fig');
