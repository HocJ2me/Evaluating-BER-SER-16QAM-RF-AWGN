function [error_rate_sig, error_rate_mut] = ofdm_awgn(SNR, M, bit_length, bit_moded, bit_sequence)

    carrier_count = 200; % S�6�3 l�0�6�6�1ng subcarrier
    symbol_count = 100; % S�6�3 l�0�6�6�1ng k�0�5 t�6�5 (symbols)
    ifft_length = 512; % �0�3�6�1 d��i c�6�5a IFFT
    CP_length = 128; % �0�3�6�1 d��i c�6�5a cyclic prefix
    CS_length = 20; % �0�3�6�1 d��i c�6�5a cyclic suffix
    alpha = 1.5/32; % H�6�3 s�6�3 c�6�5a c�6�1a s�6�7 RCOS
    carrier_position = 29:228; % V�6�7 tr�� c�6�5a subcarrier
    conj_position = 485:-1:286; % V�6�7 tr�� c�6�5a subcarrier li��n conj

    %% IFFT
    % Chuy�6�9n �0�4�6�7i t�6�9ng k�0�5 t�6�5 th��nh d�5�5ng t��n hi�6�3u th�6�5i gian
    ifft_position = zeros(ifft_length, symbol_count);
    bit_moded = reshape(bit_moded, carrier_count, symbol_count);
    ifft_position(carrier_position, :) = bit_moded(:,:);
    ifft_position(conj_position, :) = conj(bit_moded(:,:));
    signal_time = ifft(ifft_position, ifft_length);

    % Th��m cyclic prefix v�� cyclic suffix
    signal_time_C = [signal_time(end-CP_length+1:end, :); signal_time];
    signal_time_C = [signal_time_C; signal_time_C(1:CS_length, :)];

    % �0�9p d�6�3ng c�6�1a s�6�7 RCOS
    signal_window = signal_time_C .* repmat(rcoswindow(alpha, size(signal_time_C, 1)), 1, symbol_count);

    %% Truy�6�7n t��n hi�6�3u qua k��nh c�� nhi�6�1u v�� �0�4a �0�4�0�6�6�5ng
    signal_Tx = reshape(signal_window, 1, []); % Chuy�6�9n th��nh d�5�5ng th�6�5i gian �0�4�6�9 truy�6�7n
    mult_path_am = [1 0.2 0.1]; % Amplitudes c�6�5a c��c �0�4a �0�4�0�6�6�5ng
    mutt_path_time = [0 20 50]; % Th�6�5i gian tr�6�1 c�6�5a c��c �0�4a �0�4�0�6�6�5ng
    path2 = 0.2 * [zeros(1, 20) signal_Tx(1:end-20)];
    path3 = 0.1 * [zeros(1, 50) signal_Tx(1:end-50)];
    signal_Tx_mult = signal_Tx + path2 + path3; % T��n hi�6�3u sau khi tr�5�7i qua �0�4a �0�4�0�6�6�5ng

    %% Th��m nhi�6�1u AWGN
    signal_power_sig = var(signal_Tx); % C�0�0ng su�5�9t c�6�5a t��n hi�6�3u kh�0�0ng c�� �0�4a �0�4�0�6�6�5ng
    signal_power_mut = var(signal_Tx_mult); % C�0�0ng su�5�9t c�6�5a t��n hi�6�3u sau �0�4a �0�4�0�6�6�5ng
    SNR_linear = 10^(SNR/10);
    noise_power_mut = signal_power_mut / SNR_linear;
    noise_power_sig = signal_power_sig / SNR_linear;
    noise_sig = randn(size(signal_Tx)) * sqrt(noise_power_sig);
    noise_mut = randn(size(signal_Tx_mult)) * sqrt(noise_power_mut);
    Rx_data_sig = signal_Tx + noise_sig;
    Rx_data_mut = signal_Tx_mult + noise_mut;

    %% Chuy�6�9n t�6�9 d�5�5ng th�6�5i gian v�6�7 t��n hi�6�3u trong t�6�1n s�6�3
    Rx_data_mut = reshape(Rx_data_mut, ifft_length + CS_length + CP_length, []);
    Rx_data_sig = reshape(Rx_data_sig, ifft_length + CS_length + CP_length, []);

    %% Lo�5�5i b�6�1 cyclic prefix v�� cyclic suffix
    Rx_data_sig(1:CP_length, :) = [];
    Rx_data_sig(end-CS_length+1:end, :) = [];
    Rx_data_mut(1:CP_length, :) = [];
    Rx_data_mut(end-CS_length+1:end, :) = [];

    %% FFT
    fft_sig = fft(Rx_data_sig);
    fft_mut = fft(Rx_data_mut);

% Zero Forcing Filter
H = fft(eye(ifft_length + CP_length + CS_length));
H = H(carrier_position, :).';

% Apply Zero Forcing Filter directly without bsxfun
Rx_data_sig = ifft(fft(Rx_data_sig) ./ H, [], 1);
Rx_data_mut = ifft(fft(Rx_data_mut) ./ H, [], 1);



    %% L�5�9y m�6�5u t��n hi�6�3u �6�7 v�6�7 tr�� c�6�5a subcarrier
    data_sig = fft_sig(carrier_position, :);
    data_mut = fft_mut(carrier_position, :);

    %% Gi�5�7i modulasi QAM �0�4�6�9 l�5�9y l�5�5i bit
    bit_demod_sig = reshape(qamdemod(data_sig, M, 'OutputType', 'bit'), [], 1);
    bit_demod_mut = reshape(qamdemod(data_mut, M, 'OutputType', 'bit'), [], 1);

    %% T��nh t�6�1 l�6�3 l�6�9i bit
    error_bit_sig = sum(bit_demod_sig ~= bit_sequence);
    error_bit_mut = sum(bit_demod_mut ~= bit_sequence);
    error_rate_sig = error_bit_sig / bit_length;
    error_rate_mut = error_bit_mut / bit_length;
end
