results()

function [ber_data_point_at_snr, ser_data_point_at_snr] = monte_carlo(sigma, sequence_length)
    K_min = 10^5; % minimum number of errors
    K_e = 0; % number of errors
    L = 10^5; % maximum length of sequence
    k = 0; % bit counter
    no_error = 0;
    symbol_error = 0;
    
    % sequence length defined here!
    M = 2^sequence_length;
    % define our constellation
    [constellation_points, bit_sequences] = create_alphabet(sequence_length);
    seq_and_constel = associate_hamming_distance(M, constellation_points, bit_sequences);


    while K_e < K_min && k < L
        k = k + 1;
        % this is the transmitter block
        % generate a random sequence to be transmitted of length
        % sequence_length
        generated_sequence = randi([0,1], 1, sequence_length);
%         if x == 0
%             x = -1;
%         end
        % this is the modulator block where a generated sequence will be
        % modulated
        transmitted_symbol = modulator(generated_sequence, seq_and_constel);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is the channel block
        n_real = normrnd(0, sigma);
        n_imag = normrnd(0, sigma);
        received_symbol = {transmitted_symbol{1} + n_real, transmitted_symbol{2} + n_imag}; % add noise to input
        % this is the detector block
        tx_seq_hat = maximum_likelihood_demodulator(seq_and_constel, received_symbol);
        errors_i = hamming_distance_finder(tx_seq_hat, generated_sequence);
        % where errors_i is the errors made for this sequence
        if errors_i == 0
            % no error
            no_error = no_error + 1;
        else
            % error
            K_e = K_e + errors_i;
            symbol_error = symbol_error + 1;
            disp("Error")
            disp("Number of errors:")
            disp(K_e)
        end
    end
    number_of_bits_generated = k*sequence_length;
    disp("Number of times ran:")
    disp(k)
    ser_data_point_at_snr = {10*log10(1/sigma^2), symbol_error/k}; % data point in the form of 1/_sigma^2, prob of error for that sigma
    ber_data_point_at_snr = {10*log10(1/sigma^2), K_e/number_of_bits_generated};

end

function output_symbol = modulator(generated_sequence, seq_and_constel)
    % sequence needs to be in the form of i.e. [1, 0, 1]
    i = 1;
    while i <= length(seq_and_constel)
        % find said sequence in the sequence and constellation cell array
        if seq_and_constel{i}{1} == generated_sequence
            output_symbol = {seq_and_constel{i}{2}, seq_and_constel{i}{3}};
        end
        i = i + 1;
    end
end

function [constellation_points, sequences] = create_alphabet(sequence_length)
    % we must do some gray encoding
    possible_sequences = {};
    iteration = 1;
    % iterate across all possible sequences
    sequences = ff2n(sequence_length);


    % if sequence length is 3 then we have 8 symbols resulting in 8-PSK
    % lets just start with BPSK and generalize
    % 2pi(m-1)/M
    % m = 1
    % generalize for bpsk
    m = 0;
    M = 2^(sequence_length);
    constellation_points = {};
    energy = 1;
    while m < M
        phase = (2*pi*m)/M;
        constellation_points{end + 1} = {energy*cos(phase); energy*sin(phase)};
        m = m + 1;
    end
end

function hamming_distance = hamming_distance_finder(sequence1, sequence2)
    % finds the hamming distance between two different sequences
    i = 1;
    hamming_distance = 0;
    while i <= length(sequence1)
        if sequence1(i) ~= sequence2(i)
            hamming_distance = hamming_distance + 1;
        end
        i = i + 1;
    end
end

function bool_value = check_if_in(number, sequence)
% equivalent to the in function in python where we see if a number is in a
% list
    i = 1;
    bool_value = false;
    while i <= length(sequence)
        if number == sequence(i)
            bool_value = true;
            break
        end
        i = i + 1;
    end
end

function seq_and_constel = associate_hamming_distance(M, constellation_points, bit_sequences)
    % mapping is a mapping between the bit sequences and their respective
    % constellation points, this will be a cell array with 1 by 3 vector in
    % each cell, the first index is the respective sequence, while the
    % others are the corresponding real and imaginary components

    % assign a bit sequence to the first index of the constellation points
    % in other words, 000, 1, 0.
    seq_and_constel = {{bit_sequences(1, 1:end), constellation_points{1}{1}, constellation_points{1}{2}}};
    number_of_seq_to_be_mapped = M;
    sequence_indexer = 2;
    off_limits_sequences = [];
    index_constellation_points = 2;
    % keep running until there are no sequences left
    while number_of_seq_to_be_mapped ~= 0 
        
        % when we assign a sequence add that index to the
        % off_limits_sequences
        if ~check_if_in(sequence_indexer, off_limits_sequences) & hamming_distance_finder(seq_and_constel{end}{1}, bit_sequences(sequence_indexer, 1:end)) == 1
            off_limits_sequences = [off_limits_sequences, sequence_indexer];
            seq_and_constel{end + 1} = {bit_sequences(sequence_indexer, 1:end), constellation_points{index_constellation_points}{1}, constellation_points{index_constellation_points}{2}};
            index_constellation_points = index_constellation_points + 1;
        end


        % we do not want to visit the same sequence
        if sequence_indexer == M
            sequence_indexer = 2;
        else
           sequence_indexer = sequence_indexer + 1;
        end

        if length(seq_and_constel) == M
            break
        end
    end


end

function d = euclidean_distance(point1, point2)
    % where point 1 and point 2 are of the form
    % {x, y} or {real, imaginary}
    % point 2 uses 2 and 3 to ignore the sequence representation in the
    % first index
    d = sqrt((point1{1} - point2{2})^2 + (point1{2} - point2{3})^2);
end

function closest_sequence = maximum_likelihood_demodulator(seq_and_constel, received_symbol)
    i = 1;
    % go through all the possible sequences and calculate the minimum
    % distance
    % recieved symbol is of the form {real + noise, imaginary + noise
    minimum_distance = 100000;
    closest_sequence = 0;
    while i <= length(seq_and_constel)
        dist_val = euclidean_distance(received_symbol, seq_and_constel{i});
        if minimum_distance > euclidean_distance(received_symbol, seq_and_constel{i})
            closest_sequence_cons_set = seq_and_constel{i};
            minimum_distance = dist_val;
        end
        i = i + 1;
    end
    closest_sequence = closest_sequence_cons_set{1};
%     disp("Closest Sequence, and corresponding symbol real and imaginary components:")
%     disp(closest_sequence_cons_set{1})
%     disp(closest_sequence_cons_set(2:end))

end


function results()
    sigma_set = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, .9, .8, .75, .6, .5, .4, .3, .25, .15, .1, .05, .01};
    
    sequence_length = {1, 2, 3};
    sequence_length_index = 1;
    ber_fxn_b = [];
    ser_fxn_b = [];
    sigma_index = 1;
    while sigma_index <= length(sigma_set)
        [ber_data_point_at_snr, ser_data_point_at_snr] = monte_carlo(sigma_set{sigma_index}, 1) % returns one data point at a time for each snr
        ber_fxn_b = [ber_fxn_b, ber_data_point_at_snr]
        ser_fxn_b = [ser_fxn_b, ser_data_point_at_snr]
        sigma_index = sigma_index + 1;    
    end

    ber_fxn_q = [];
    ser_fxn_q = [];
    sigma_index = 1;
    while sigma_index <= length(sigma_set)
        [ber_data_point_at_snr, ser_data_point_at_snr] = monte_carlo(sigma_set{sigma_index}, 2) % returns one data point at a time for each snr
        ber_fxn_q = [ber_fxn_q, ber_data_point_at_snr]
        ser_fxn_q = [ser_fxn_q, ser_data_point_at_snr]
        sigma_index = sigma_index + 1;    
    end

    ber_fxn_8 = [];
    ser_fxn_8 = [];
    sigma_index = 1;
    while sigma_index <= length(sigma_set)
        [ber_data_point_at_snr, ser_data_point_at_snr] = monte_carlo(sigma_set{sigma_index}, 3) % returns one data point at a time for each snr
        ber_fxn_8 = [ber_fxn_8, ber_data_point_at_snr]
        ser_fxn_8 = [ser_fxn_8, ser_data_point_at_snr]
        sigma_index = sigma_index + 1;    
    end

    ber_x_b = []
    ber_y_b = []
    ser_x_b = []
    ser_y_b = []
    i = 1;
    while i <= length(ber_fxn_b) - 1
        ber_x_b = [ber_x_b, ber_fxn_b{i}]
        ber_y_b = [ber_y_b, ber_fxn_b{i+1}]
        ser_x_b = [ser_x_b, ser_fxn_b{i}]
        ser_y_b = [ser_y_b, ser_fxn_b{i+1}]
        i = i + 2
    end

    ber_x_q = []
    ber_y_q = []
    ser_x_q = []
    ser_y_q = []
    i = 1;
    while i <= length(ber_fxn_q) - 1
        ber_x_q = [ber_x_q, ber_fxn_q{i}]
        ber_y_q = [ber_y_q, ber_fxn_q{i+1}]
        ser_x_q = [ser_x_q, ser_fxn_q{i}]
        ser_y_q = [ser_y_q, ser_fxn_q{i+1}]
        i = i + 2
    end

    ber_x_8 = []
    ber_y_8 = []
    ser_x_8 = []
    ser_y_8 = []
    i = 1;
    while i <= length(ber_fxn_8) - 1
        ber_x_8 = [ber_x_8, ber_fxn_8{i}]
        ber_y_8 = [ber_y_8, ber_fxn_8{i+1}]
        ser_x_8 = [ser_x_8, ser_fxn_8{i}]
        ser_y_8 = [ser_y_8, ser_fxn_8{i+1}]
        i = i + 2
    end

            display(ber_x)
            display(ber_y)
    semilogy(ber_x_b,ber_y_b,'LineWidth',2.0)
    hold on
    semilogy(ber_x_q,ber_y_q,'LineWidth',2.0)
    hold on
    semilogy(ber_x_8,ber_y_8,'LineWidth',2.0)
    hold on
    xlim([-15 30])
    ylim([10^-10 1])
    legend('BPSK Simulated', 'QPSK Simulated', '8PSK Simulated')
    xlabel('SNR (dB)')
    ylabel('BER')
    set(gcf,'color','w');
    set(gca,'Color','w');
    title("Simulated Bit Error Rate")
    grid
    graph the most recent data
    sequence_length_index = sequence_length_index + 1;
    disp(sigma_index)

    hold off
    semilogy(ser_x_b, ser_y_b,'LineWidth',2.0)
    hold on
    semilogy(ser_x_q, ser_y_q,'LineWidth',2.0)
    hold on
    semilogy(ser_x_8, ser_y_8,'LineWidth',2.0)
    hold on
    xlim([-15 30])
    ylim([10^-10 1])
    legend('BPSK Simulated', 'QPSK Simulated', '8PSK Simulated')
    xlabel('SNR (dB)')
    ylabel('SER')
    set(gcf,'color','w');
    set(gca,'Color','w');
    graph the most recent data
    sequence_length_index = sequence_length_index + 1;
    title("Simulated Symbol Error Rate")
    grid

    EbNo = (-15:25)';
    [berB, serB] = berawgn(EbNo,'psk',2,'nondiff');
    [berQ, serQ] = berawgn(EbNo,'psk',4,'nondiff');
    [ber8, ser8] = berawgn(EbNo,'psk',8,'nondiff');
    semilogy(EbNo, berB, 'LineWidth', 1.5)
    hold on
    semilogy(EbNo, [berQ ber8], 'LineWidth', 1.0)
    hold off
    xlim([-15 20])
    ylim([10^-10 1])
    xlabel('SNR (dB)')
    ylabel('BER')
    legend('BPSK Theoretical','QPSK Theoretical','8PSK Theoretical')
    set(gcf,'color','w');
    set(gca,'Color','w');
    title("Theoretical Bit Error Rate")
    grid
    EbNo = (-15:25)';
    [berB, serB] = berawgn(EbNo,'psk',2,'nondiff');
    [berQ, serQ] = berawgn(EbNo,'psk',4,'nondiff');
    [ber8, ser8] = berawgn(EbNo,'psk',8,'nondiff');
    semilogy(EbNo, serB, 'LineWidth', 1.5)
    hold on
    semilogy(EbNo, [serQ ser8], 'LineWidth', 1.5)
    hold off
    xlim([-15 20])
    ylim([10^-10 1])
    xlabel('SNR (dB)')
    ylabel('SER')
    set(gcf,'color','w');
    set(gca,'Color','w');
    legend('BPSK Theoretical','QPSK Theoretical','8PSK Theoretical')
    title("Theoretical Symbol Error Rate")
    grid
end