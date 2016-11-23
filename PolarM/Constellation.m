classdef Constellation < handle
    % Deals with constellation transmit and receive
    
    properties
        
        constellation_points;
        
        n_bits;
        n_sym;
        bit_sym_map;

        constellation_name;
        is_gray;
        
    end
    
    properties (Constant)
        
        bpsk = [1 -1];
        
        ask4_gray = [-3, -1, 3, 1]/sqrt(5);
        
        ask4_sp = (-3:2:3)/sqrt(5);

        ask8_gray = [-7, -5, -1, -3, 7, 5, 1, 3]/sqrt(21);
        
        ask8_sp = (-7:2:7)/sqrt(21);
        
        ask16_gray = [-15, -13, -9, -11, -1, -3, -7, -5, ...
                        15, 13, 9, 11, 1, 3, 7, 5]/sqrt(85);

       ask16_sp = (-15:2:15)/sqrt(85);
        
    end
    
    methods
        
        function obj = Constellation ( constellion_name)
            obj.constellation_name = constellion_name;
            obj.is_gray = 0;
            if strcmp (constellion_name, 'bpsk')
                obj.constellation_points = obj.bpsk;
                obj.n_bits = 1;
            elseif strcmp (constellion_name, 'ask4-gray')
                obj.constellation_points = obj.ask4_gray;
                obj.n_bits = 2;
                obj.is_gray = 1;
            elseif strcmp (constellion_name, 'ask4-sp')
                obj.constellation_points = obj.ask4_sp;
                obj.n_bits = 2;
            elseif strcmp (constellion_name, 'ask8-gray')
                obj.constellation_points = obj.ask8_gray;
                obj.n_bits = 3;
                obj.is_gray = 1;
            elseif strcmp (constellion_name, 'ask8-sp')
                obj.constellation_points = obj.ask8_sp;
                obj.n_bits = 3;
            elseif strcmp (constellion_name, 'ask16-gray')
                obj.constellation_points = obj.ask16_gray;
                obj.n_bits = 4;
                obj.is_gray = 1;               
            elseif strcmp (constellion_name, 'ask16-sp')
                obj.constellation_points = obj.ask16_sp;
                obj.n_bits = 4;
            else
                disp('Unsupported constellation');
            end
            
            obj.n_sym = 2^obj.n_bits;
            obj.bit_sym_map = zeros(obj.n_sym, obj.n_bits);
            for sym_index = 0 : obj.n_sym - 1
                a = dec2bin(sym_index, obj.n_bits);
                for bit_index = 1 :obj.n_bits
                    if a(obj.n_bits + 1 - bit_index) == '1'
                        obj.bit_sym_map(sym_index + 1, bit_index) = 1;
                    end
                end
            end
            
            obj.constellation_points = obj.constellation_points/sqrt(mean(obj.constellation_points.^2));
                        
        end
        
        function [sym] = modulate(obj, bits)
            sym = zeros(floor(length(bits)/obj.n_bits), 1);
            for i = 1 : 1 : length(bits)/obj.n_bits
                symbol = 0;
                for j = 1 : obj.n_bits
                    symbol = symbol + 2^(j-1) * bits((i-1)*obj.n_bits+j);
                end
                sym(i) = obj.constellation_points(symbol + 1);                    
            end
        end
        
        function [p1, llr] = compute_llr_mlc(obj, y, n_0, u)
            p0 = zeros(length(y), 1);
            p1 = zeros(length(y), 1);
            for y_index = 1 : length(y)
                for sym_index = 1  : obj.n_sym
                    valid_symbol = 1;
                    for m_index = 1 : size(u, 2)
                        if obj.bit_sym_map(sym_index, m_index) ~= u(y_index, m_index)
                            valid_symbol = 0;
                        end
                    end
                    if valid_symbol == 0
                        continue;
                    end
                    m_index = size(u, 2) + 1;
                    if obj.bit_sym_map(sym_index, m_index) == 0
                        p0(y_index) = p0(y_index) + exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0);
                    else
                        p1(y_index) = p1(y_index) + exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0);
                    end
                    
                end
            end
            llr = log(p0./p1);
            p1 = p1./(p0 + p1);
            
        end

        function [p1, llr] = compute_llr_bicm(obj, y, n_0)
            p0 = zeros(length(y) * obj.n_bits, 1);
            p1 = zeros(length(y) * obj.n_bits, 1);
            for y_index = 1 : length(y)
                for sym_index = 1  : 2^obj.n_bits
                    if length(n_0) == 1
                        p_sym = exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0);
                    else
                        p_sym = exp(-abs(y(y_index) - obj.constellation_points(sym_index))^2/2/n_0(y_index));
                    end
                    for m_index = 1 : obj.n_bits
                        if obj.bit_sym_map(sym_index, m_index) == 0
                            p0((y_index-1) * obj.n_bits +  m_index) = p0((y_index-1) * obj.n_bits +  m_index) + p_sym;
                        else
                            p1((y_index-1) * obj.n_bits +  m_index) = p1((y_index-1) * obj.n_bits +  m_index) + p_sym;
                        end
                    end
                end
            end
            llr = log(p0./p1);
            p1 = p1./(p0 + p1);
        end
                       
        function mean_llr = get_mean_llr_bicm(obj, snr_db)
            num_sym = 10000;
            bits = rand(obj.n_bits*num_sym, 1) < 0.5;
            mod_sym = obj.modulate(bits);
            
            noise = randn(length(mod_sym), 1);
            sigma = sqrt(1/2) *  10^(-snr_db/20);
            y = mod_sym + noise * sigma;

            [~, llr_bicm] =  obj.compute_p1_bicm(y, sigma^2);
            
            llr_bicm = llr_bicm .* (1 - 2*bits);
            mean_llr = zeros(obj.n_bits, 1);
            for i_bit = 1 : obj.n_bits
                mean_llr(i_bit) = mean(llr_bicm(i_bit: obj.n_bits:end));
            end
        end

        function mean_llr = get_mean_llr_mlc(obj, snr_db)
            num_sym = 20000;
            bits = rand(obj.n_bits*num_sym, 1) < 0.5;
            mod_sym = obj.modulate(bits);
            noise = randn(length(mod_sym), 1);

            sigma = sqrt(1/2) *  10^(-snr_db/20);
            y = mod_sym + noise * sigma;
            
            llr_mlc = zeros(num_sym*obj.n_bits, 1);
            u = zeros(num_sym, obj.n_bits);
            n_0 = sigma^2;
            for i_bit = 1 : obj.n_bits
                u(:, i_bit) = bits(i_bit:obj.n_bits:end);
                [~, llr_mlc(i_bit:obj.n_bits:end)] = obj.compute_p1_mlc(y, n_0, u(:, 1:i_bit-1));
            end
            
            llr_mlc = llr_mlc .* (1 - 2*bits);
            mean_llr = zeros(obj.n_bits, 1);
            for i_bit = 1 : obj.n_bits
                mean_llr(i_bit) = mean(llr_mlc(i_bit: obj.n_bits:end));
            end
            
        end

        
        function cap_vec = get_mlc_capacity(obj, snr_db)
            sigma = sqrt(1/2) *  10^(-snr_db/20);
            n_0 = sigma^2;
            
            y_min = max(obj.constellation_points) + 6*sigma + 1;
            delta_y = sigma * 0.01;
            y_vec = -y_min: delta_y: y_min;
            
            h_y = zeros(obj.n_bits, 1);
            h_y_given_u = zeros(obj.n_bits, 1);
            
            for i_bit = 1 : obj.n_bits
                num_sets = 2^(i_bit - 1);
                
                p_y = zeros(length(y_vec), num_sets);
                p_y_u = zeros(length(y_vec), num_sets, 2);
                
                for y_index = 1 : length(y_vec)
                    y = y_vec(y_index);
                    for i_sym = 1 : obj.n_sym
                        x = obj.bit_sym_map(i_sym);
                        set_index = 1;
                        if i_bit > 1
                            a = dec2bin(i_sym - 1, obj.n_bits);
                            a = a(obj.n_bits - i_bit + 2:obj.n_bits);
                            set_index = bin2dec(a) + 1;
                        end
                        p_y(y_index, set_index) = p_y(y_index, set_index) + ...
                            exp(-(y - obj.constellation_points(i_sym))^2/2/n_0)/sqrt(2*pi*n_0)/obj.n_sym ;
                        p_y_u(y_index, set_index, obj.bit_sym_map(i_sym, i_bit) + 1 ) = ...
                            p_y_u(y_index, set_index, obj.bit_sym_map(i_sym, i_bit) + 1) + ...
                            exp(-(y - obj.constellation_points(i_sym))^2/2/n_0)/sqrt(2*pi*n_0)/obj.n_sym * 2;
                    end
                end
                for y_index = 1 : length(y_vec)
                    for set_index = 1 : num_sets
                        if p_y(y_index, set_index) > 0
                            h_y(i_bit) = h_y(i_bit) + (-log2(p_y(y_index, set_index)))*p_y(y_index, set_index) * delta_y;
                        end
                    end
                end
                
                for y_index = 1 : length(y_vec)
                    for set_index = 1 : num_sets
                        for u = 1 : 2
                            if p_y_u(y_index, set_index, u) > 0
                                h_y_given_u(i_bit) = h_y_given_u(i_bit) + ...
                                    (-log2(p_y_u(y_index, set_index, u))) * ...
                                    p_y_u(y_index,set_index, u) * 0.5  * delta_y;
                            end
                        end
                    end
                end
                
            end
            cap_vec = h_y - h_y_given_u;
            
        end
        
        
        function cap_vec = get_bicm_capacity(obj, snr_db)

            sigma = sqrt(1/2) *  10^(-snr_db/20);
            n_0 = sigma^2;
            y_min = max(obj.constellation_points) + 6*sigma + 1;
            delta_y = sigma * 0.1;
            y_vec = -y_min: delta_y: y_min;
            cap_vec = zeros(obj.n_bits, 1);
            for i_bit = 1 : obj.n_bits
                p_y = zeros(length(y_vec), 1);
                p_y_u = zeros(length(y_vec), 2);
                for y_index = 1 : length(y_vec)
                    y = y_vec(y_index);
                    for i_sym = 1 : obj.n_sym
                        p_y(y_index) = p_y(y_index) + ...
                            exp(-(y - obj.constellation_points(i_sym))^2/2/n_0)/sqrt(2*pi*n_0)/obj.n_sym ;
                        p_y_u(y_index, obj.bit_sym_map(i_sym, i_bit) + 1 ) = ...
                            p_y_u(y_index, obj.bit_sym_map(i_sym, i_bit) + 1) + ...
                            exp(-(y - obj.constellation_points(i_sym))^2/2/n_0)/sqrt(2*pi*n_0)/obj.n_sym * 2;
                    end
                end               
                h_y = 0;
                h_y_u = 0;
                for y_index = 1 : length(y_vec)
                    if p_y(y_index) > 0
                        h_y     = h_y + log2(p_y(y_index))*p_y(y_index)*delta_y * (-1);
                    end
                    if p_y_u(y_index, 1) > 0
                        h_y_u     = h_y_u + 0.5 * log2(p_y_u(y_index, 1))*p_y_u(y_index, 1)*delta_y * (-1);
                    end
                    if p_y_u(y_index, 2) > 0
                        h_y_u     = h_y_u + 0.5 * log2(p_y_u(y_index, 2))*p_y_u(y_index, 2)*delta_y * (-1);
                    end
                end
                cap_vec(i_bit) = h_y - h_y_u;
            end
        end
        
        function cap_vec = get_polarized_capacity(obj, snr_db)
            
            mat_file_name = ['CapacityHelper/PolarizedCapacityData/', obj.constellation_name, '_snr_', num2str(snr_db), '.mat'];
            if exist(mat_file_name, 'file')
               cap_data = load(mat_file_name);
               cap_vec = cap_data.cap_vec;
            else
                disp(['Determining and storing polarized capacity for constellation = ' , obj.constellation_name, ' for snr = ', num2str(snr_db)]);
                num_sym = 250e3;
                num_codes = obj.n_bits;
                if num_codes == 3
                    num_codes = 4;
                end
                info_bits = rand(obj.n_bits*num_sym, 1) < 0.5;
                coded_bits = zeros(obj.n_bits*num_sym, 1);
                
                for i_sym = 1 : floor(obj.n_bits*num_sym/num_codes)
                    coded_bits((i_sym-1)*num_codes + 1 : i_sym * num_codes) = ...
                        PolarCode.polar_encode(info_bits((i_sym-1)*num_codes + 1 : i_sym * num_codes));
                end
                
                mod_sym = obj.modulate(coded_bits);
                
                noise = randn(length(mod_sym), 1);
                sigma = sqrt(1/2) *  10^(-snr_db/20);
                y = mod_sym + noise * sigma;
                
                [~, llr_bicm] =  obj.compute_llr_bicm(y, sigma^2);
                
                u_llr = zeros(obj.n_bits*num_sym, 1);
                
                for i_sym = 1 : floor(obj.n_bits*num_sym/num_codes)
                    [~, u_llr((i_sym-1)*num_codes + 1 : i_sym * num_codes)] = ...
                        PolarCode.polar_decode_capacity_llr(llr_bicm((i_sym-1)*num_codes + 1 : i_sym * num_codes), ...
                        info_bits((i_sym-1)*num_codes + 1 : i_sym * num_codes));
                end
                
                y_max = 100;
                delta_y = 0.25;
                y_vec = -y_max:delta_y:y_max;
                
                cap_vec = zeros(num_codes, 1);
                u_llr = max(u_llr, -y_max);
                u_llr = min(u_llr, y_max);
                
                for i_bit = 1 : num_codes
                    p_y = zeros(length(y_vec), 1);
                    p_y_u = zeros(length(y_vec), 2);
                    num_sym0 = 0;
                    num_sym1 = 0;
                    for i_sym = 1 : floor(obj.n_bits*num_sym/num_codes)
                        y_index = floor( (u_llr(i_bit + (i_sym-1)*num_codes) + y_max )/delta_y) + 1;
                        p_y(y_index) = p_y(y_index) + 1;
                        if info_bits(i_bit + (i_sym-1)*num_codes) == 0
                            p_y_u(y_index, 1) = p_y_u(y_index, 1) + 1;
                            num_sym0 = num_sym0 + 1;
                        else
                            p_y_u(y_index, 2) = p_y_u(y_index, 2) + 1;
                            num_sym1 = num_sym1 + 1;
                            
                        end
                    end
                    p_y = p_y /floor(obj.n_bits*num_sym/num_codes);
                    p_y_u(:,1) = p_y_u(:,1)/num_sym0;
                    p_y_u(:,2) = p_y_u(:,2)/num_sym1;
                    h_y = 0;
                    h_y_u = 0;
                    for y_index = 1 : length(y_vec)
                        if p_y(y_index) > 0
                            h_y     = h_y + log2(p_y(y_index))*p_y(y_index) * (-1);
                        end
                        if p_y_u(y_index, 1) > 0
                            h_y_u     = h_y_u + 0.5 * log2(p_y_u(y_index, 1))*p_y_u(y_index, 1) * (-1);
                        end
                        if p_y_u(y_index, 2) > 0
                            h_y_u     = h_y_u + 0.5 * log2(p_y_u(y_index, 2))*p_y_u(y_index, 2) * (-1);
                        end
                    end
                    cap_vec(i_bit) = h_y - h_y_u;
                    
                end
                cap_vec = min(cap_vec, 1);
                save(mat_file_name, 'cap_vec');
            end
        end
    end       
end

