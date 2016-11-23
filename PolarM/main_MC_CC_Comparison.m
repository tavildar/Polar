

%main_MC_CC_Comparison

% close all;
% clear all;
n = 10;
block_length = 2^n;
info_length = block_length/2;

constellation_vec = {'ask4-sp', 'ask4-gray', 'ask8-gray', 'ask16-sp', 'ask16-gray'};

snrMap = containers.Map(...
    {'bpsk', 'ask4-sp', 'ask4-gray', 'ask8-gray', 'ask16-sp', 'ask16-gray'}, ...
    {-0.5, 4.5, 4.5, 9, 12, 13});

rxMap = containers.Map(...
    {'bpsk', 'ask4-sp', 'ask4-gray', 'ask8-gray', 'ask16-sp', 'ask16-gray'}, ...
    {'bicm', 'mlc' , 'bicm', 'bicm', 'mlc', 'bicm'});

monte_carlo_runs = 250e3;

polar_code = PolarCode(block_length, info_length, 0.5);

max_runs = 100e3;
max_err = 250;

base_snr_db_vec = -3: 0.25: 3;
% num_block_err = zeros(length(base_snr_db_vec), length(constellation_vec));
% num_runs = zeros(length(base_snr_db_vec), length(constellation_vec));

figure(1); 
tic
for i_const = 3 : length(constellation_vec)
    const_name = constellation_vec{i_const};
    receiver_algo = rxMap(const_name);
    design_snr_db = snrMap(const_name);
    modulation = Constellation(const_name);
    
    polar_code.monte_carlo_code_construction(design_snr_db, monte_carlo_runs, const_name, receiver_algo);
    
    snr_db_vec = design_snr_db  + base_snr_db_vec;
    
    for i_run = 1 : max_runs
        
        if mod(i_run, ceil(max_runs/10)) == 1
            disp(['Sim iteration running = ', num2str(i_run)]);
        end
        info = rand(1 , info_length) > 0.5;
        
        if strcmp(receiver_algo, 'bicm')
            coded_bits = polar_code.encode(info);
        elseif strcmp(receiver_algo, 'mlc')
            num_codes = modulation.n_bits;
            info_paded = zeros(block_length, 1);
            coded_bits = zeros(block_length, 1);
            info_paded(polar_code.info_bits) = info;
            for layer = 1 : num_codes
                info_bit_loc = (layer-1) * block_length/num_codes + 1 : layer * block_length/num_codes;
                coded_bits(layer:num_codes:block_length) = PolarCode.polar_encode(info_paded(info_bit_loc))';
            end
        end
        
        mod_sym = modulation.modulate(coded_bits);
        
        prev_decoded = zeros(length(snr_db_vec), 1);
        
        noise = randn(length(mod_sym), 1);
        
        
        for i_snr = 1 : length(snr_db_vec)
            
            if num_block_err(i_snr, i_const) > max_err
                continue;
            end
            run_sim = 1;
            num_runs(i_snr, i_const) =  num_runs(i_snr, i_const) + 1;
            for i_snr2 = 1 : i_snr
                if prev_decoded(i_snr2)
                    run_sim = 0;
                end
            end
            
            if (run_sim == 0)
                continue;
            end
            
            snr_db = snr_db_vec(i_snr);
            
            sigma = sqrt(1/2) *  10^(-snr_db/20);
            
            y = mod_sym + noise * sigma;
            if strcmp(receiver_algo, 'bicm')
                p1 = 0.5*ones(block_length, 1);
                p1(1:floor(block_length/modulation.n_bits) * modulation.n_bits) =  modulation.compute_llr_bicm(y, sigma^2);
                decoded_bits = polar_code.decode_sc_p1(p1');
            else
                decoded_info = zeros(block_length/num_codes, num_codes);
                decoded_coded = zeros(block_length/num_codes, num_codes);
                for layer = 1 : num_codes
                    u1 = decoded_coded(:, 1:layer-1);
                    p1 = modulation.compute_llr_mlc(y, sigma^2, u1).';
                    info_bit_loc = (layer-1) * block_length/num_codes + 1 : layer * block_length/num_codes;
                    [decoded_info(:, layer), decoded_coded(:, layer)] = PolarCode.polar_decode(p1, polar_code.frozen_bits(info_bit_loc));
                end
                u_decoded = decoded_info(:);
                decoded_bits = u_decoded(polar_code.info_bits)';
            end
            
            err = any(info ~= decoded_bits);
            
            if err
                num_block_err(i_snr, i_const) =  num_block_err(i_snr, i_const) + 1;
            else
                prev_decoded(i_snr) = 1;
            end
            
        end
    end
    
    semilogy(snr_db_vec + 10*log10(block_length/info_length) - 10*log10(modulation.n_bits), num_block_err(:, i_const)./num_runs(:,i_const), 'LineWidth' , 2);
    hold on; grid on;
end
toc
legend(constellation_vec);
savefig('monte-carlo.fig');
