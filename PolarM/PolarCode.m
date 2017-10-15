%
% Created by Saurabh Tavildar on 3/23/16.
%
classdef PolarCode < handle
    % Polar code
    properties
        
        block_length;
        info_length;
        frozen_bits;
        n;
        design_epsilon;
        
        bit_reversed_order;
        info_bits;
        
        %sc_decoder
        u_sc;
        x_sc;
        
        %scl decoder
        list_size;
        p_scl;
        c_scl;
        i_scl;
        
        inactivePathIndices;
        inactivePathIndicesSize;
        activePathArray;
        pathIndexToArrayIndex;
        inactiveArrayIndices;
        inactiveArrayIndicesSize;
        arrayReferenceCount;
        
        lambda_offset;
        list_offset;
        
        llr_based_computation;
        llr_scl;
        llr_path_metric;
        
        % CRC
        crc_matrix;
        crc_size;
        
        %code construction
        cc_method;
        cc_parameter;
        cc_misc;
        
        %
        info_bit_order
    end
    
    methods
        
        
        %% code construction
        function obj = PolarCode(block_length, info_length, design_epsilon, crc_size)
            
            obj.block_length = block_length;
            obj.info_length = info_length;
            obj.n = log2(block_length);
            obj.design_epsilon = design_epsilon;
            
            if nargin < 4
                crc_size = 0;
            end
            
            obj.crc_size = crc_size;
            
            if obj.crc_size ~= 0
                obj.crc_matrix = floor(2*rand(obj.crc_size, info_length));
            end
            
            obj.bit_reversed_order = bitrevorder((1:obj.block_length)');
            
            channels = PolarCode.calculate_channel_polarization(design_epsilon, obj.n );
            channels = channels(obj.bit_reversed_order);
            [~, info_bits_sorted] = sort(channels, 'ascend');
            obj.info_bit_order = info_bits_sorted(obj.block_length:-1:1);
            obj.frozen_bits = ones(1,obj.block_length);
            obj.frozen_bits(info_bits_sorted(1:obj.info_length + obj.crc_size)) = 0;
            obj.info_bits = info_bits_sorted(1:obj.info_length + obj.crc_size);
            obj.llr_based_computation = 0;
            
            obj.cc_method = 'bhattacharya';
            obj.cc_parameter = design_epsilon;
            obj.cc_misc = '';
            
            disp(['Bhattacharya code construction done. BLER estimate:', num2str(sum(channels(obj.info_bits)))]);
            
        end
        
        function monte_carlo_code_construction(obj, design_snr_db, num_runs, constellation_name, receiver_algo)
            if (nargin < 3) || isempty(num_runs)
                num_runs = 100e3;
            end
            if (nargin < 4) || isempty(constellation_name)
                constellation_name = 'bpsk';
            end
            if (nargin < 5) || isempty(receiver_algo)
                receiver_algo = 'bicm';
            end
            
            
            obj.cc_method = 'monte-carlo';
            obj.cc_parameter = design_snr_db;
            obj.cc_misc = [constellation_name, '_', receiver_algo, '_', num2str(num_runs)];
            
            txt_file_name = ['CodeConstructionData/MC_block_length_', obj.get_unique_string(), '.txt'];
            
            if exist(txt_file_name, 'file')
                channels = load(txt_file_name);
                disp(['Found monte carlo BER pattern: ', txt_file_name]);
            else
                disp(['Didnt find monte carlo BER here: ', txt_file_name]);
                disp('Hence, determining the  monte carlo ber.');
                channels = obj.monte_carlo(design_snr_db, num_runs, constellation_name, receiver_algo);
                fileID = fopen(txt_file_name,'w');
                for c  = channels
                    fprintf(fileID,'%d \n',c);
                end
                fclose(fileID);
            end
            [~, channel_order] = sort(channels, 'ascend');
            if channels(channel_order(obj.info_length + obj.crc_size)) < 100
                disp('Warning: not enough runs to get a reliable code.');
                disp(['Worst # of errors = ', num2str(channels(channel_order(obj.info_length + obj.crc_size)))]);
            end
            
            obj.info_bits = channel_order(1:obj.info_length + obj.crc_size);
            
            obj.frozen_bits = ones(1,obj.block_length);
            obj.frozen_bits(obj.info_bits) = 0;
            bler_estimate = sum(channels(obj.info_bits))/num_runs;
            
            
            disp(['Monte carlo code construction done. Bler estimate = ', num2str(bler_estimate)]);
            
        end
        
        function num_err = monte_carlo(obj, design_snr_db, num_runs, constellation_name, receiver)
            
            num_err = zeros(obj.block_length, 1);
            modulation = Constellation(constellation_name);
            
            for i_run  = 1 : num_runs
                if mod(i_run, ceil(num_runs/10)) == 1
                    disp(['Code construction iteration running = ', num2str(i_run)]);
                end
                if strcmp(receiver, 'bicm')
                    
                    dummy_info = (rand(1, obj.block_length) < 0.5);
                    dummy_coded = PolarCode.polar_encode(dummy_info);
                    
                elseif strcmp(receiver, 'mlc')
                    
                    num_codes = modulation.n_bits;
                    dummy_info = (rand(num_codes, obj.block_length/num_codes) < 0.5);
                    dummy_coded = zeros(obj.block_length, 1);
                    for layer = 1 : num_codes
                        dummy_coded(layer:num_codes:obj.block_length) = PolarCode.polar_encode(dummy_info(layer,:))';
                    end
                    
                end
                
                mod_sym = modulation.modulate(dummy_coded);
                sigma = sqrt(1/2) *  10^(-design_snr_db/20);
                noise = sigma * randn(length(mod_sym), 1);
                y = mod_sym + noise;
                
                if strcmp(receiver, 'bicm')
                    
                    p1 = 0.5 * ones(obj.block_length, 1);
                    eff_block_length = floor(obj.block_length/modulation.n_bits)*modulation.n_bits;
                    [p1(1:eff_block_length),~] = modulation.compute_llr_bicm(y, sigma^2);
                    [~, ber_tmp] = PolarCode.polar_decode_monte(p1, dummy_info);
                    num_err = num_err + ber_tmp';
                    
                elseif strcmp(receiver, 'mlc')
                    
                    decoded_coded = zeros(obj.block_length/num_codes, num_codes);
                    ber = zeros(obj.block_length/num_codes, num_codes);
                    
                    for layer = 1 : num_codes
                        u = decoded_coded(:, 1:layer-1);
                        [p1,~] = modulation.compute_llr_mlc(y, sigma^2, u);
                        [decoded_coded(:, layer), ber(:, layer)] = PolarCode.polar_decode_monte(p1, dummy_info(layer, :));
                    end
                    ber = ber(:);
                    num_err = num_err + ber;
                    
                end
            end
        end
        
        function [bler_estimate] = ga_code_construction(obj, design_snr_db, constellation_name, receiver_algo)
            
            modulation = Constellation(constellation_name);
            
            if strcmp(receiver_algo, 'bicm')
                
                if modulation.n_bits == 1
                    capacity = modulation.get_bicm_capacity(design_snr_db);
                else
                    capacity = modulation.get_polarized_capacity(design_snr_db);
                end
                
            elseif strcmp(receiver_algo, 'mlc')
                
                capacity = modulation.get_mlc_capacity(design_snr_db);
                
            end
            
            num_codes = modulation.n_bits;
            
            mean_llrs = get_bpsk_llr_for_capacity(capacity);
            
            llr_vec = zeros(obj.block_length, 1);
            
            if modulation.n_bits == 3  && strcmp(receiver_algo, 'mlc') %non power of 2 not supported
                disp('MLC code construction not supported');
            else
                for i_code = 1 : num_codes
                    bit_loc = (i_code-1) * obj.block_length/num_codes + 1 : i_code * obj.block_length/num_codes;
                    llr_vec(bit_loc) = mean_llrs(i_code);
                end
            end
            
            
            bit_rev_order_subcode = bitrevorder(1:obj.block_length/ num_codes);
            channels = zeros(obj.block_length, 1);
            
            for i_code = 1 : num_codes
                bit_loc = (i_code-1) * obj.block_length/num_codes + 1 : i_code * obj.block_length/num_codes;
                tmp_c = calculate_awgn_polarization(llr_vec((bit_loc)).', obj.n - log2(num_codes));
                channels(bit_loc) = tmp_c(bit_rev_order_subcode);
            end
            
            [~, channel_order] = sort(channels, 'descend');
            obj.info_bit_order = channel_order(obj.block_length:-1:1);
            
            obj.info_bits = channel_order(1:obj.info_length + obj.crc_size);
            obj.frozen_bits = ones(1,obj.block_length);
            obj.frozen_bits(obj.info_bits) = 0;
            bler_estimate =  sum(qfunc(  sqrt(channels(obj.info_bits))/sqrt(2)));
            
            obj.cc_method = 'gauss-approx';
            obj.cc_parameter = design_snr_db;
            obj.cc_misc = [constellation_name, '_', receiver_algo];
            
            disp(['GA code construction done. BLER estimate:', num2str(bler_estimate)]);
            
        end
        
        
        function [unique_string] = get_unique_string(obj)
            unique_string = [num2str(obj.block_length), '_', num2str(length(obj.info_bits)),  ...
                '_cc_method_', obj.cc_method, '_cc_param_', num2str(obj.cc_parameter), '_', obj.cc_misc];
        end
        
        
        %% encoder
        
        function coded_bits = encode(obj, info_bits)
            
            info_bits_padded = zeros(obj.block_length, 1);
            if obj.crc_size ~= 0
                crc = mod(obj.crc_matrix*info_bits', 2);
                info_bits = [info_bits, crc'];
            end
            info_bits_padded(obj.info_bits) = info_bits;
            coded_bits = PolarCode.polar_encode(info_bits_padded);
            
        end
        
        
        %%  Decoder helper
        
        function [b] = crc_check(obj, info_bits)
            info = info_bits(1:obj.info_length);
            crc  = info_bits(obj.info_length + 1:obj.info_length + obj.crc_size);
            crc_check =  mod(obj.crc_matrix*info', 2)';
            b = 1 - any(crc ~= crc_check);
        end
        
        %%  SC decoder
        
        function decoded_bits = decode_sc_p1(obj, p1)
            
            [decoded_bits, ~] = PolarCode.polar_decode(p1, obj.frozen_bits);
            decoded_bits = decoded_bits(obj.info_bits(1:obj.info_length));
            
        end
        
        %%  SCL decoder
        
        function [u] = decode_scl_p1(obj, p1, p0, list_size)
            
            obj.list_size = list_size;
            obj.llr_based_computation  = 0;
            obj.initializeDataStructures();
            l_index = obj.assignInitialPath();
            s_index = obj.getArrayPointer_P(0, l_index);
            obj.p_scl(obj.get_i_scl(0, 0, s_index) + 1:obj.get_i_scl(0, obj.block_length - 1, s_index) + 1 , 1) = p0;
            obj.p_scl(obj.get_i_scl(0, 0, s_index) + 1:obj.get_i_scl(0, obj.block_length - 1, s_index) + 1 , 2) = p1;
            u = obj.polar_decode_scl();
            
        end
        
        function [u] = decode_scl_llr(obj, llr, list_size)
            
            obj.list_size = list_size;
            obj.llr_based_computation = 1;
            obj.initializeDataStructures();
            l_index = obj.assignInitialPath();
            s_index = obj.getArrayPointer_P(0, l_index);
            obj.llr_scl( obj.get_i_scl(0, 0, s_index) + 1:obj.get_i_scl(0, obj.block_length - 1, s_index) + 1) = llr;
            u = obj.polar_decode_scl();
            
        end
        
        function [u] = polar_decode_scl(obj)
            
            for phi = 0 : obj.block_length - 1
                
                obj.recursivelyCalcP_scl(obj.n, phi);
                
                if obj.frozen_bits(phi + 1) == 1
                    obj.continuePaths_FrozenBit(phi);
                else
                    obj.continuePaths_UnfrozenBit(phi);
                end
                
                if mod(phi, 2) == 1
                    obj.recursivelyUpdateC_scl(obj.n, phi);
                end
                
            end
            
            l_index = obj.findMostProbablePath(1);
            c_m = obj.getArrayPointer_C(obj.n, l_index);
            info = obj.i_scl(c_m+1,:);
            u = info(obj.info_bits(1:obj.info_length));
            
        end
        
        function index = get_i_scl(obj, lambda, beta, list_index)
            index = beta + obj.lambda_offset(lambda + 1) + obj.list_offset(list_index + 1);
        end
        
        function init_i_scl(obj)
            obj.lambda_offset = (2.^(obj.n - (0: obj.n)) - 1);
            obj.list_offset = (0:obj.list_size)*(2 * obj.block_length - 1);
            
        end
        
        function initializeDataStructures(obj)
            
            obj.inactivePathIndices = zeros(obj.list_size,1);
            obj.inactivePathIndicesSize = 0;
            % the above two variables are used to define a stack
            
            obj.activePathArray =  zeros(obj.list_size,1);
            obj.pathIndexToArrayIndex = zeros(obj.n  + 1, obj.list_size);
            
            obj.inactiveArrayIndices = zeros(obj.n  + 1, obj.list_size);
            obj.inactiveArrayIndicesSize = zeros(obj.n + 1, 1);
            % the above two variables are used to define a vector of stacks
            
            obj.arrayReferenceCount = zeros(obj.n  + 1, obj.list_size);
            
            if obj.llr_based_computation
                obj.llr_scl = zeros(obj.list_size * (2 * obj.block_length - 1), 1);
                obj.llr_path_metric =  zeros(obj.list_size, 1);
            else
                obj.p_scl = zeros(obj.list_size * (2 * obj.block_length - 1), 2);
            end
            
            obj.c_scl = zeros(obj.list_size * (2 * obj.block_length - 1), 2);
            obj.i_scl = zeros(obj.list_size, obj.block_length);
            obj.init_i_scl();
            
            for lambda = 0 : obj.n
                for i_list = 0 : obj.list_size - 1
                    obj.inactiveArrayIndices(lambda + 1, i_list + 1) = i_list;
                    
                end
                obj.inactiveArrayIndicesSize(lambda + 1) = obj.list_size;
            end
            
            for i_list = 0 : obj.list_size - 1
                obj.activePathArray(i_list + 1) = 0;
                obj.inactivePathIndices(i_list + 1) = i_list;
            end
            
            obj.inactivePathIndicesSize  = obj.list_size;
            
        end
        
        function l_index = assignInitialPath(obj)
            l_index = obj.inactivePathIndices(obj.inactivePathIndicesSize);
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize - 1;
            obj.activePathArray(l_index + 1) = 1;
            
            for lambda = 0 : obj.n
                s = obj.inactiveArrayIndices(lambda + 1, obj.inactiveArrayIndicesSize(lambda + 1));
                obj.inactiveArrayIndicesSize(lambda + 1) = obj.inactiveArrayIndicesSize(lambda + 1) - 1;
                obj.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s;
                obj.arrayReferenceCount(lambda + 1, l_index + 1) = 1;
            end
        end
        
        function l_p_index = clonePath(obj, l_index)
            l_p_index = obj.inactivePathIndices(obj.inactivePathIndicesSize);
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize - 1;
            obj.activePathArray(l_p_index + 1) = 1;
            
            if obj.llr_based_computation
                obj.llr_path_metric(l_p_index + 1) = obj.llr_path_metric(l_index + 1);
            end
            
            for lambda = 0 : obj.n
                s = obj.pathIndexToArrayIndex(lambda + 1, l_index + 1);
                obj.pathIndexToArrayIndex(lambda + 1, l_p_index + 1) = s;
                obj.arrayReferenceCount(lambda + 1, s + 1) = obj.arrayReferenceCount(lambda + 1, s + 1) + 1;
            end
        end
        
        function killPath(obj, l_index)
            obj.activePathArray(l_index + 1) = 0;
            obj.inactivePathIndices(obj.inactivePathIndicesSize + 1) = l_index;
            obj.inactivePathIndicesSize = obj.inactivePathIndicesSize  + 1;
            if obj.llr_based_computation
                obj.llr_path_metric(l_index + 1) = 0;
            end
            for lambda = 0 : obj.n
                s = obj.pathIndexToArrayIndex(lambda + 1, l_index + 1);
                obj.arrayReferenceCount(lambda + 1, s + 1) = obj.arrayReferenceCount(lambda + 1, s + 1) - 1;
                if obj.arrayReferenceCount(lambda + 1, s + 1) == 0
                    obj.inactiveArrayIndices(lambda + 1, obj.inactiveArrayIndicesSize(lambda + 1) + 1) = s;
                    obj.inactiveArrayIndicesSize(lambda + 1)  = obj.inactiveArrayIndicesSize(lambda + 1)  + 1;
                end
            end
        end
        
        function [b] = pathIndexInactive(obj, l_index)
            b = 1 - obj.activePathArray(l_index + 1);
        end
        
        function [s_p] = getArrayPointer_P(obj, lambda, l_index)
            s = obj.pathIndexToArrayIndex(lambda + 1, l_index + 1);
            m = obj.n;
            if obj.arrayReferenceCount(lambda + 1, s + 1) == 1
                s_p = s;
            else
                s_p = obj.inactiveArrayIndices(lambda + 1, obj.inactiveArrayIndicesSize(lambda + 1));
                i_s_p = obj.lambda_offset(lambda + 1) + obj.list_offset(s_p + 1) + 1: ...
                    obj.lambda_offset(lambda + 1) + obj.list_offset(s_p + 1) + 2^(m - lambda);
                i_s = obj.lambda_offset(lambda + 1) + obj.list_offset(s + 1) + 1: ...
                    obj.lambda_offset(lambda + 1) + obj.list_offset(s + 1) + 2^(m - lambda);
                obj.c_scl(i_s_p, :) = obj.c_scl(i_s, :);
                if obj.llr_based_computation
                    obj.llr_scl(i_s_p) = obj.llr_scl(i_s);
                else
                    obj.p_scl(i_s_p, :) = obj.p_scl(i_s, :);
                end
                obj.inactiveArrayIndicesSize(lambda + 1) = obj.inactiveArrayIndicesSize(lambda + 1) - 1;
                obj.arrayReferenceCount(lambda + 1, s + 1) = obj.arrayReferenceCount(lambda + 1, s + 1) - 1;
                obj.arrayReferenceCount(lambda + 1, s_p + 1) = 1;
                obj.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s_p;
            end
        end
        
        function [s_p] = getArrayPointer_C(obj, lambda, l_index)
            s = obj.pathIndexToArrayIndex(lambda + 1, l_index + 1);
            m = obj.n;
            if obj.arrayReferenceCount(lambda + 1, s + 1) == 1
                s_p = s;
            else
                s_p = obj.inactiveArrayIndices(lambda + 1, obj.inactiveArrayIndicesSize(lambda + 1));
                i_s_p = obj.lambda_offset(lambda + 1) + obj.list_offset(s_p + 1) + 1: ...
                    obj.lambda_offset(lambda + 1) + obj.list_offset(s_p + 1) + 2^(m - lambda);
                i_s = obj.lambda_offset(lambda + 1) + obj.list_offset(s + 1) + 1: ...
                    obj.lambda_offset(lambda + 1) + obj.list_offset(s + 1) + 2^(m - lambda);
                if obj.llr_based_computation
                    obj.llr_scl(i_s_p) = obj.llr_scl(i_s);
                else
                    obj.p_scl(i_s_p, :) = obj.p_scl(i_s, :);
                end
                obj.c_scl(i_s_p, :) = obj.c_scl(i_s, :);
                obj.inactiveArrayIndicesSize(lambda + 1) = obj.inactiveArrayIndicesSize(lambda + 1) - 1;
                obj.arrayReferenceCount(lambda + 1, s + 1) = obj.arrayReferenceCount(lambda + 1, s + 1) - 1;
                obj.arrayReferenceCount(lambda + 1, s_p + 1) = 1;
                obj.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s_p;
            end
        end
        
        
        function recursivelyCalcP_scl(obj, lambda, phi)
            
            if lambda == 0
                return;
            end
            
            psi = floor(phi/2);
            if mod(phi, 2) == 0
                obj.recursivelyCalcP_scl(lambda - 1, psi);
            end
            
            sigma  = 0;
            p_index_3_base_list = zeros(obj.list_size, 1);
            l_index_1_list = zeros(obj.list_size, 1);
            for l_index = 0 : obj.list_size - 1
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                l_index_1_list(l_index+1) = obj.getArrayPointer_P(lambda, l_index);
                l_index_2 = obj.getArrayPointer_P(lambda - 1, l_index);
                l_index_3 = obj.getArrayPointer_C(lambda, l_index);
                
                p_index_1_base = obj.lambda_offset(lambda) + obj.list_offset(l_index_2 + 1) + 1;
                p_index_3_base_list(l_index+1) = obj.lambda_offset(lambda + 1) + obj.list_offset(l_index_1_list(l_index+1) + 1) + 1;
                c_index_3_base = obj.lambda_offset(lambda + 1) + obj.list_offset(l_index_3 + 1) + 1;
                for beta = 0: 2^(obj.n - lambda) - 1
                    p_index_1 = p_index_1_base + 2 * beta;
                    p_index_2 = p_index_1_base + 2 * beta + 1;
                    p_index_3 = p_index_3_base_list(l_index+1) + beta;
                    if mod(phi, 2) == 0
                        if obj.llr_based_computation
                            if max( abs(obj.llr_scl ( p_index_1)), abs(obj.llr_scl ( p_index_2)) ) < 40
                                obj.llr_scl(p_index_3) = PolarCode.cnop_llr( obj.llr_scl ( p_index_1), obj.llr_scl ( p_index_2));
                                % log( (exp( obj.llr_scl ( p_index_1) + obj.llr_scl ( p_index_2)) + 1) ...
                                %    /(exp( obj.llr_scl ( p_index_1))  + exp( obj.llr_scl ( p_index_2))) );
                            else
                                obj.llr_scl(p_index_3) = sign( obj.llr_scl ( p_index_1)) * sign(obj.llr_scl ( p_index_2)) * min(abs(obj.llr_scl ( p_index_2)), abs(obj.llr_scl ( p_index_1)));
                            end
                        else
                            obj.p_scl(p_index_3, 1) =  0.5 * ( obj.p_scl ( p_index_1 , 1) *  obj.p_scl(p_index_2, 1)  +  obj.p_scl ( p_index_1 , 2) *  obj.p_scl(p_index_2, 2));
                            obj.p_scl(p_index_3, 2) =  0.5 * ( obj.p_scl ( p_index_1 , 2) *  obj.p_scl(p_index_2, 1)  +  obj.p_scl ( p_index_1 , 1) *  obj.p_scl(p_index_2, 2));
                            sigma = max(obj.p_scl(p_index_3, 1), sigma);
                            sigma = max(obj.p_scl(p_index_3, 2), sigma);
                        end
                    else
                        u_p = obj.c_scl( c_index_3_base + beta, 1);
                        if obj.llr_based_computation
                            obj.llr_scl(p_index_3) = (-1)^u_p * obj.llr_scl(p_index_1) +  obj.llr_scl(p_index_2);
                        else
                            obj.p_scl(p_index_3, 1)  = 0.5 * obj.p_scl (p_index_1, mod(u_p, 2) + 1) * ...
                                obj.p_scl(p_index_2, 1);
                            obj.p_scl(p_index_3, 2)  = 0.5 * obj.p_scl (p_index_1, mod(u_p + 1, 2) + 1) * ...
                                obj.p_scl(p_index_2, 2);
                            sigma = max(obj.p_scl(p_index_3, 1), sigma);
                            sigma = max(obj.p_scl(p_index_3, 2), sigma);
                        end
                    end
                end
            end
            
            for l_index = 0 : obj.list_size - 1
                if sigma == 0 %typically happens because of underflow
                    break;
                end
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                p_range = p_index_3_base_list(l_index+1): p_index_3_base_list(l_index+1) +  2^(obj.n - lambda) - 1;
                obj.p_scl(p_range, :) = obj.p_scl(p_range, :)/sigma;
            end
            
        end
        
        
        function recursivelyUpdateC_scl(obj, lambda, phi)
            if mod(phi, 2) == 0
                disp('Error: phi should always be odd in this function call');
            end
            psi = floor(phi/2);
            
            for l_index = 0 : obj.list_size - 1
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                l_index_1 = obj.getArrayPointer_C(lambda, l_index);
                l_index_2 = obj.getArrayPointer_C(lambda - 1, l_index);
                
                p_index_1 = obj.lambda_offset(lambda + 1) + obj.list_offset(l_index_1 + 1)  + 1;
                p_index_2 = obj.lambda_offset(lambda) + obj.list_offset(l_index_2 + 1)  +  1;
                
                for beta = 0: 2^(obj.n - lambda) - 1
                    obj.c_scl(p_index_2 + 2*beta, mod(psi, 2) + 1) = ...
                        mod( obj.c_scl(p_index_1 + beta, 1) + ...
                        obj.c_scl(p_index_1 + beta, 2),  2);
                    obj.c_scl(p_index_2 + 2 * beta + 1, mod(psi, 2) + 1) = ...
                        obj.c_scl(p_index_1 + beta, 2);
                end
                
            end
            
            if mod(psi, 2) == 1
                obj.recursivelyUpdateC_scl(lambda - 1, psi);
            end
            
        end
        
        function continuePaths_FrozenBit(obj, phi)
            
            for l_index = 0 : obj.list_size -1
                
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                l_index_1 = obj.getArrayPointer_C(obj.n, l_index);
                obj.c_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, mod(phi, 2) + 1) = 0;
                if obj.llr_based_computation
                    obj.llr_path_metric(l_index_1 + 1) = obj.llr_path_metric(l_index_1 + 1) ...
                        + log(1 + exp(-obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1)));
                end
            end
            
        end
        
        function    continuePaths_UnfrozenBit(obj, phi)
            
            probForks = -realmax * ones(obj.list_size, 2);
            index = 0;
            for l_index = 0 : obj.list_size - 1
                
                if obj.activePathArray(l_index + 1)
                    l_index_1 = obj.getArrayPointer_P(obj.n, l_index);
                    if obj.llr_based_computation
                        % computing negative of path metric so that an
                        % ascending order can be used for sorting
                        probForks(l_index + 1, 1) =  - (obj.llr_path_metric(l_index_1 + 1) ...
                            + log(1 + exp(-obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1))));
                        probForks(l_index + 1, 2) =  - ( obj.llr_path_metric(l_index_1 + 1) ...
                            + log(1 + exp(obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1))));
                    else
                        probForks(l_index + 1, 1) = obj.p_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, 1);
                        probForks(l_index + 1, 2) = obj.p_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, 2);
                    end
                    index = index + 1;
                    
                end
            end
            
            rho = min(2*index, obj.list_size);
            contForks = zeros(obj.list_size, 2);
            prob = sort(probForks(:), 'descend');
            
            threshold = prob(rho);
            num_populated = 0;
            for l_index = 0 : obj.list_size - 1
                for j_index = 1 : 2
                    if num_populated == rho
                        break;
                    end
                    if  probForks(l_index + 1, j_index) > threshold
                        contForks(l_index + 1, j_index) = 1;
                        num_populated = num_populated + 1;
                    end
                end
            end
            
            if num_populated < rho
                for l_index = 0 : obj.list_size - 1
                    for j_index = 1 : 2
                        if num_populated == rho
                            break;
                        end
                        if  probForks(l_index + 1, j_index) == threshold
                            contForks(l_index + 1, j_index) = 1;
                            num_populated = num_populated + 1;
                        end
                    end
                end
            end
            
            
            for l_index = 0 : obj.list_size - 1
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                if (contForks(l_index + 1, 1) == 0) && ( contForks(l_index + 1, 2) == 0)
                    obj.killPath(l_index);
                end
            end
            
            for l_index = 0 : obj.list_size - 1
                if (contForks(l_index + 1, 1) == 0) && ( contForks(l_index + 1, 2) == 0)
                    continue;
                end
                l_index_1 = obj.getArrayPointer_C(obj.n, l_index);
                if (contForks(l_index + 1, 1) == 1) && ( contForks(l_index + 1, 2) == 1)
                    obj.c_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 0;
                    obj.i_scl(l_index_1 + 1, phi + 1) = 0;
                    
                    l_p = obj.clonePath(l_index);
                    
                    l_index_2 = obj.getArrayPointer_C(obj.n, l_p);
                    obj.i_scl(l_index_2 + 1, 1 : phi) = obj.i_scl(l_index_1 + 1, 1 : phi);
                    obj.c_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_2 + 1) + 1, mod(phi,2) + 1) = 1;
                    obj.i_scl(l_index_2 + 1, phi + 1) = 1;
                    if obj.llr_based_computation
                        obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) ...
                            + log(1 + exp(-obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1)));
                        obj.llr_path_metric(l_p + 1) = obj.llr_path_metric(l_p + 1) ...
                            + log(1 + exp(obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_2 + 1) + 1)));
                    end
                    
                else
                    if contForks(l_index + 1, 1) == 1
                        obj.c_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 0;
                        obj.i_scl(l_index_1 + 1, phi + 1) = 0;
                        if obj.llr_based_computation
                            obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) ...
                                + log(1 + exp(-obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1)));
                        end
                    else
                        obj.c_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 1;
                        obj.i_scl(l_index_1 + 1, phi + 1) = 1;
                        if obj.llr_based_computation
                            obj.llr_path_metric(l_index + 1) = obj.llr_path_metric(l_index + 1) ...
                                + log(1 + exp(obj.llr_scl(obj.lambda_offset(obj.n + 1) + obj.list_offset(l_index_1 + 1) + 1)));
                        end
                    end
                end
                
            end
            
        end
        
        function [l_p_index] = findMostProbablePath(obj, crc_check)
            l_p_index = 0;
            
            p_max = 0;
            if obj.llr_based_computation
                p_max = realmax;
            end
            path_with_crc = 0;
            for l_index = 0 : obj.list_size -1
                
                if obj.activePathArray(l_index + 1) == 0
                    continue;
                end
                
                c_index = obj.getArrayPointer_C( obj.n, l_index);
                if (crc_check) && (obj.crc_size ~= 0)
                    a = obj.i_scl(c_index+1,:);
                    u = a(obj.info_bits);
                    if obj.crc_check(u) == 0
                        continue;
                    end
                end
                path_with_crc = 1;
                if obj.llr_based_computation
                    if p_max > obj.llr_path_metric(l_index + 1)
                        p_max = obj.llr_path_metric(l_index + 1);
                        l_p_index = l_index;
                    end
                else
                    p_index = obj.getArrayPointer_P( obj.n, l_index);
                    if p_max < obj.p_scl(obj.get_i_scl( obj.n, 0, p_index) + 1, obj.c_scl(obj.get_i_scl( obj.n, 0, c_index) + 1, 2) + 1)
                        l_p_index = l_index;
                        p_max  = obj.p_scl(obj.get_i_scl( obj.n, 0, p_index) + 1, obj.c_scl(obj.get_i_scl( obj.n, 0, c_index) + 1, 2) + 1);
                    end
                end
            end
            
            if (crc_check) && (path_with_crc == 0) % no path with crc check found
                l_p_index = obj.findMostProbablePath(0);
            end
            
        end
        
        %% helper function
        function [ bler, ber ] = get_bler_quick(obj, ebno_vec, list_size_vec)
            
            snr_db_vec = ebno_vec + 10*log10(obj.info_length/obj.block_length);
            
            num_block_err = zeros(length(ebno_vec), length(list_size_vec));
            num_bit_err = zeros(length(ebno_vec), length(list_size_vec));
            num_runs = zeros(length(ebno_vec), length(list_size_vec));
            max_err = 50;
            max_runs = 500;
            
            for i_run = 1 : max_runs
                
                if mod(i_run, ceil(max_runs/10)) == 1
                    disp(['Sim iteration running = ', num2str(i_run)]);
                end
                info = rand(1 , obj.info_length) < 0.5;
                coded_bits = obj.encode(info);
                bpsk = 2 * coded_bits - 1;
                sigma = sqrt(1/2);
                noise = sigma * randn(1, obj.block_length);
                prev_decoded = zeros(length(list_size_vec), length(ebno_vec));
                
                for i_ebno = 1 : length(ebno_vec)
                    snr_db = snr_db_vec(i_ebno);
                    received_bits = 10^(snr_db/20) * bpsk + noise;
                    p1 = exp(-(received_bits - 10^(snr_db/20)).^2/(2 * sigma^2))/sigma/sqrt(2*pi);
                    p0 = exp(-(received_bits + 10^(snr_db/20)).^2/(2 * sigma^2))/sigma/sqrt(2*pi);
                    
                    for i_list = 1 : length(list_size_vec)
                        if num_block_err(i_ebno, i_list) > max_err
                            continue;
                        end
                        
                        num_runs(i_ebno,  i_list) =  num_runs(i_ebno, i_list) + 1;
                        
                        run_sim = 1;
                        
                        for i_ebno2 = 1 : i_ebno
                            if prev_decoded(i_list, i_ebno2)
                                run_sim = 0;
                            end
                        end
                        
                        if (run_sim == 0)
                            % This is a hack to speed up simulations --
                            % it assumes that this run will be decoded correctly since it was
                            % decoded correctly for a lower EbNo
                            continue;
                        end
                        if list_size_vec(i_list) == 1
                            decoded_bits = obj.decode_sc_p1(p1./(p1+p0));
                        else
                            %  decoded_bits = obj.decode_scl_p1(p1, p0, list_size_vec(i_list));
                            decoded_bits = obj.decode_scl_llr(log(p0./p1), list_size_vec(i_list));
                        end
                        err = any(info ~= decoded_bits);
                        if err
                            num_block_err(i_ebno, i_list) =  num_block_err(i_ebno, i_list) + 1;
                            num_bit_err(i_ebno, i_list) = num_bit_err(i_ebno, i_list) + sum(info ~= decoded_bits);
                        else
                            prev_decoded(i_list, i_ebno) = 1;
                        end
                    end
                end
                
            end
            bler = num_block_err./num_runs;
            ber = num_bit_err./num_runs;
            
        end
        
    end
    
    methods (Static)
        % the next four functions are modified version of code in:
        % http://pfister.ee.duke.edu/courses/ecen655/polar.pdf
        function x = polar_encode(info_bits_padded)
            
            N = length(info_bits_padded);
            if (N == 1)
                x = info_bits_padded;
            else
                u1u2 = mod(info_bits_padded(1:2:end) + info_bits_padded(2:2:end) , 2);
                u2 = info_bits_padded(2:2:end);
                x = [PolarCode.polar_encode(u1u2) PolarCode.polar_encode(u2)];
            end
            
        end
        
        function [u, x] = polar_decode(y,f)
            N = length(y);
            if (N==1)
                if (f == 0)  % info bit
                    x = (1-sign(1-2*y))/2;
                else  % frozen bit -- assumed to be zero
                    x = 0;
                end
                u = x;
            else
                u1est = PolarCode.cnop(y(1:2:end),y(2:2:end));
                [uhat1,u1hardprev] = PolarCode.polar_decode(u1est,f(1:N/2));
                u2est = PolarCode.vnop(PolarCode.cnop(u1hardprev,y(1:2:end)),y(2:2:end));
                [uhat2,u2hardprev] = PolarCode.polar_decode(u2est,f(N/2+1:end));
                u = [uhat1 uhat2];
                x = reshape([PolarCode.cnop(u1hardprev,u2hardprev); u2hardprev],1,[]);
            end
        end
        
        function z = cnop(w1,w2) % notation = probability of bit = 1
            z = w1.*(1-w2) + w2.*(1-w1);
        end
        
        function z = vnop(w1,w2) % notation = probability of bit = 1
            z = w1.*w2 ./ (w1.*w2 + (1-w1).*(1-w2));
        end
        
        function [x, ber] = polar_decode_monte(y, dummy_info)
            N = length(y);
            if (N==1)
                if (y > 0.5 && dummy_info == 1) || (y <= 0.5 && dummy_info == 0)
                    ber = 0;
                else
                    ber = 1;
                end
                x = dummy_info;
            else
                u1est = PolarCode.cnop(y(1:2:end),y(2:2:end));
                [u1hardprev, ber1] = PolarCode.polar_decode_monte(u1est,dummy_info(1:N/2));
                u2est = PolarCode.vnop(PolarCode.cnop(u1hardprev',y(1:2:end)),y(2:2:end));
                [u2hardprev, ber2] = PolarCode.polar_decode_monte(u2est,dummy_info(N/2+1:end));
                ber = [ber1 ber2];
                x = reshape([PolarCode.cnop(u1hardprev,u2hardprev); u2hardprev],1,[]);
            end
        end
        
        function [x, u] = polar_decode_capacity(y, dummy_info)
            N = length(y);
            if (N==1)
                u = y;
                x = dummy_info;
            else
                u1est = PolarCode.cnop(y(1:2:end),y(2:2:end));
                [u1hardprev, u1] = PolarCode.polar_decode_capacity(u1est,dummy_info(1:N/2));
                u2est = PolarCode.vnop(PolarCode.cnop(u1hardprev',y(1:2:end)),y(2:2:end));
                [u2hardprev, u2] = PolarCode.polar_decode_capacity(u2est,dummy_info(N/2+1:end));
                u = [u1 u2];
                x = reshape([PolarCode.cnop(u1hardprev,u2hardprev); u2hardprev],1,[]);
            end
        end
        
        
        function [x, u] = polar_decode_capacity_llr(y, dummy_info)
            N = length(y);
            if (N==1)
                u = y;
                x = dummy_info;
            else
                u1est = PolarCode.cnop_llr(y(1:2:end),y(2:2:end));
                [u1hardprev, u1] = PolarCode.polar_decode_capacity_llr(u1est,dummy_info(1:N/2));
                u2est = PolarCode.vnop_llr((1 - 2*u1hardprev').*y(1:2:end),y(2:2:end));
                [u2hardprev, u2] = PolarCode.polar_decode_capacity_llr(u2est,dummy_info(N/2+1:end));
                u = [u1 u2];
                x = reshape([PolarCode.cnop(u1hardprev,u2hardprev); u2hardprev],1,[]);
            end
        end
        
        
        function l = cnop_llr(l1,l2)
            l = 2 * atanh(tanh(l1/2).*tanh(l2/2));
        end
        
        function l = vnop_llr(l1,l2)
            l = l1 + l2;
        end
        
        
        function channels  = calculate_channel_polarization( epsilon, n)
            
            N  = 2^n;
            log_domain_enabled = 1;
            if log_domain_enabled
                epsilon = log(epsilon);
            end
            if length(epsilon) == 1
                channels = epsilon*ones(1, N);
            elseif length(epsilon) == N
                channels = epsilon;
            end
            
            for i = 1:n
                c1 = channels(1:2:N);
                c2 = channels(2:2:N);
                if log_domain_enabled
                    channels = [log(exp(c1) + exp(c2) - exp(c1+c2)), c1 + c2];
                else
                    channels = [c1 + c2 - c1.*c2, c1.*c2];
                end
            end
            
            if log_domain_enabled
                channels = exp(channels);
            end
            
        end
        
    end
    
end


