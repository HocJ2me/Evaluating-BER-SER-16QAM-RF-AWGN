function [error_rate_sig, error_rate_mut] = ofdm_awgn(SNR, M, bit_length, bit_moded, bit_sequence)

    carrier_count = 200; % Số lượng subcarrier
    symbol_count = 100; % Số lượng ký tự (symbols)
    ifft_length = 512; % Độ dài của IFFT
    CP_length = 128; % Độ dài của cyclic prefix
    CS_length = 20; % Độ dài của cyclic suffix
    alpha = 1.5/32; % Hệ số của cửa sổ RCOS
    carrier_position = 29:228; % Vị trí của subcarrier
    conj_position = 485:-1:286; % Vị trí của subcarrier liên conj

    %% IFFT
    % Chuyển đổi từng ký tự thành dạng tín hiệu thời gian
    ifft_position = zeros(ifft_length, symbol_count);
    bit_moded = reshape(bit_moded, carrier_count, symbol_count);
    ifft_position(carrier_position, :) = bit_moded(:,:);
    ifft_position(conj_position, :) = conj(bit_moded(:,:));
    signal_time = ifft(ifft_position, ifft_length);

    % Thêm cyclic prefix và cyclic suffix
    signal_time_C = [signal_time(end-CP_length+1:end, :); signal_time];
    signal_time_C = [signal_time_C; signal_time_C(1:CS_length, :)];

    % Áp dụng cửa sổ RCOS
    signal_window = signal_time_C .* repmat(rcoswindow(alpha, size(signal_time_C, 1)), 1, symbol_count);

    %% Truyền tín hiệu qua kênh có nhiễu và đa đường
    signal_Tx = reshape(signal_window, 1, []); % Chuyển thành dạng thời gian để truyền
    mult_path_am = [1 0.2 0.1]; % Amplitudes của các đa đường
    mutt_path_time = [0 20 50]; % Thời gian trễ của các đa đường
    path2 = 0.2 * [zeros(1, 20) signal_Tx(1:end-20)];
    path3 = 0.1 * [zeros(1, 50) signal_Tx(1:end-50)];
    signal_Tx_mult = signal_Tx + path2 + path3; % Tín hiệu sau khi trải qua đa đường

    %% Thêm nhiễu AWGN
    signal_power_sig = var(signal_Tx); % Công suất của tín hiệu không có đa đường
    signal_power_mut = var(signal_Tx_mult); % Công suất của tín hiệu sau đa đường
    SNR_linear = 10^(SNR/10);
    noise_power_mut = signal_power_mut / SNR_linear;
    noise_power_sig = signal_power_sig / SNR_linear;
    noise_sig = randn(size(signal_Tx)) * sqrt(noise_power_sig);
    noise_mut = randn(size(signal_Tx_mult)) * sqrt(noise_power_mut);
    Rx_data_sig = signal_Tx + noise_sig;
    Rx_data_mut = signal_Tx_mult + noise_mut;

    %% Chuyển từ dạng thời gian về tín hiệu trong tần số
    Rx_data_mut = reshape(Rx_data_mut, ifft_length + CS_length + CP_length, []);
    Rx_data_sig = reshape(Rx_data_sig, ifft_length + CS_length + CP_length, []);

    %% Loại bỏ cyclic prefix và cyclic suffix
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



    %% Lấy mẫu tín hiệu ở vị trí của subcarrier
    data_sig = fft_sig(carrier_position, :);
    data_mut = fft_mut(carrier_position, :);

    %% Giải modulasi QAM để lấy lại bit
    bit_demod_sig = reshape(qamdemod(data_sig, M, 'OutputType', 'bit'), [], 1);
    bit_demod_mut = reshape(qamdemod(data_mut, M, 'OutputType', 'bit'), [], 1);

    %% Tính tỷ lệ lỗi bit
    error_bit_sig = sum(bit_demod_sig ~= bit_sequence);
    error_bit_mut = sum(bit_demod_mut ~= bit_sequence);
    error_rate_sig = error_bit_sig / bit_length;
    error_rate_mut = error_bit_mut / bit_length;
end
