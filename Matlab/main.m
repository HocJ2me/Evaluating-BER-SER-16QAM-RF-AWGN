function ofdm_rayleigh_zf_16qam_simulation_plot()
    % Thiết lập tham số
    SNR_dB_range = 0:2:20; % Phạm vi giả định SNR từ 0 đến 20 dB
    num_simulations = length(SNR_dB_range);
    ber_array = zeros(1, num_simulations);
    ser_array = zeros(1, num_simulations);
    
     
    % Mô phỏng từng giá trị SNR
    for i = 1:num_simulations
        SNR_dB = SNR_dB_range(i);
        [ber, ser] = ofdm_rayleigh_zf_16qam_simulation(SNR_dB);
        ber_array(i) = ber;
        ser_array(i) = ser;
    end
    
    % Vẽ biểu đồ
    figure;
    semilogy(SNR_dB_range, ber_array, 'o-', 'LineWidth', 2, 'DisplayName', 'BER');
    hold on;
    semilogy(SNR_dB_range, ser_array, 's-', 'LineWidth', 2, 'DisplayName', 'SER');
    grid on;
    title('OFDM Simulation with Rayleigh Channel and Zero Forcing Equalization (16QAM)');
    xlabel('SNR (dB)');
    ylabel('Error Rate');
    ylim([0.4, 1.1]);
    legend('Location', 'best');
end

function [ber, ser] = ofdm_rayleigh_zf_16qam_simulation(SNR_dB)
    % Thiết lập tham số OFDM
    N = 64; % Số lượng subcarrier
    M = 16; % Số mức của modulasi (16QAM)
    num_symbols = 100; % Số lượng ký tự (symbols)
    cp_length = 16; % Độ dài của cyclic prefix
    
    % Tạo ngẫu nhiên dãy bit để modulasi
    input_bits = randi([0 1], N * num_symbols * log2(M), 1);
    
    % Modulasi 16QAM
    modulated_symbols = qammod(input_bits, M, 'InputType', 'bit');
    
    % Chuyển đổi thành ma trận ký tự (symbols)
    input_matrix = reshape(modulated_symbols, N, num_symbols);
    
    % Thực hiện IFFT
    time_domain_symbols = ifft(input_matrix);
    
    % Thêm cyclic prefix
    time_domain_symbols_with_cp = [time_domain_symbols(end - cp_length + 1:end, :); time_domain_symbols];
    
    % Truyền qua kênh Rayleigh với nhiễu trắng AWGN
    SNR_linear = 10^(SNR_dB / 10);
    channel_coefficients = (randn(1, num_symbols) + 1i * randn(1, num_symbols)) / sqrt(2); % Kênh Rayleigh
    received_symbols = time_domain_symbols_with_cp .* repmat(channel_coefficients, N + cp_length, 1);
    noise_power = 1 / SNR_linear;
    noise = (randn(N + cp_length, num_symbols) + 1i * randn(N + cp_length, num_symbols)) * sqrt(noise_power);
    received_symbols_with_noise = received_symbols + noise;
    
    % Loại bỏ cyclic prefix
    received_symbols_no_cp = received_symbols_with_noise(cp_length + 1:end, :);
    
    % Thực hiện FFT
    frequency_domain_symbols = fft(received_symbols_no_cp);

    
    % Zero Forcing Equalization
    H = fft(eye(N));
    H = H(1:N, :); % Adjust the size of H to match frequency_domain_symbols
    H = H / norm(H, 'fro'); % Normalize to avoid noise amplification
    
    equalized_symbols = zeros(size(frequency_domain_symbols));
    for col = 1:num_symbols
        H_col = H(:, mod(col - 1, size(H, 2)) + 1);
        equalized_symbols(:, col) = frequency_domain_symbols(:, col) ./ H_col;
    end

    % Giải modulasi QAM để lấy lại bit
    demodulated_bits = qamdemod(equalized_symbols(:), M, 'OutputType', 'bit');
    
    % Tính tỷ lệ lỗi bit
    ber = sum(demodulated_bits ~= input_bits) / length(input_bits);
 
    % Tính Symbol Error Rate (SER)
    symbols_original = qammod(input_bits, M, 'InputType', 'bit');
    symbols_received = qammod(demodulated_bits, M, 'InputType', 'bit');
    ser = sum(symbols_original ~= symbols_received) / length(symbols_original);

    
end
