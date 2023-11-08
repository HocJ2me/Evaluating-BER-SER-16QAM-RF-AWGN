%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share — to copy, distribute and transmit the work
% to Remix — to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating 16-QAM transmission and reception and compare the 
% simulated and theoretical symbol error probability

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 9 December 2007
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% symbol error rate for 16-QAM modulation

% symbol error rate for 16-QAM modulation
clear
N = 2*10^5; % number of symbols
alpha16qam = [-3 -1 1 3]; % 16-QAM alphabets
Es_N0_dB = [0:20]; % multiple Es/N0 values
ipHat = zeros(1,N);
for ii = 1:length(Es_N0_dB)
    ip = randsrc(1,N,alpha16qam) + j*randsrc(1,N,alpha16qam);
    s = (1/sqrt(10))*ip; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    y_re = real(y); % real part
    y_im = imag(y); % imaginary part

    ipHat_re(find(y_re< -2/sqrt(10)))           = -3;
    ipHat_re(find(y_re > 2/sqrt(10)))           =  3;
    ipHat_re(find(y_re>-2/sqrt(10) & y_re<=0))  = -1;
    ipHat_re(find(y_re>0 & y_re<=2/sqrt(10)))   =  1;

    ipHat_im(find(y_im< -2/sqrt(10)))           = -3;
    ipHat_im(find(y_im > 2/sqrt(10)))           =  3;
    ipHat_im(find(y_im>-2/sqrt(10) & y_im<=0))  = -1;
    ipHat_im(find(y_im>0 & y_im<=2/sqrt(10)))   =  1;
    ipHat = ipHat_re + j*ipHat_im;
    nErr(ii) = size(find([ip- ipHat]),2); % couting the number of errors
end

simBer = nErr/N;
theoryBer = 3/2*erfc(sqrt(0.1*(10.^(Es_N0_dB/10))));
close all
figure
semilogy(Es_N0_dB,theoryBer,'b.-','LineWidth',2);
hold on
semilogy(Es_N0_dB,simBer,'mx-','Linewidth',2);
axis([0 20 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16-QAM modulation') 












