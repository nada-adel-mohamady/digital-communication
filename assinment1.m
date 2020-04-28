pkg load communications
close all
N = 5 %number of bits to be sent 
%data generation
sent_bits = randi([0 1],1,N);

%---if sent_bit=0 make it -1 otherwise =1
S = 2*sent_bits-1;
SS = 2*sent_bits-1;
A = 1  % amplitude 
T = 1  % time interval , sample at T
%---x axis ------
x_axis = linspace(0,1,100);
g=ones(1,T)* A; %rectangular pulse
%-------------------------------
fs = 100;
S = upsample(S,fs);
h1 = ones(1,fs);
tmp = conv(S,h1);
sqrFilter = ones(1,fs)./fs;
triFilter = linspace(0,1,fs)./fs;
noFilter = [1 zeros(1,fs-1)];
%---------------------------------
%---------PLOTING----------------
z1 = conv(tmp,sqrFilter);
z1 = downsample(z1,fs);
z1 = z1(2:end);
z2 = conv(tmp,noFilter);
z3 = conv(tmp,triFilter);
figure('Name','Case1 Matched Filter','NumberTitle','off')
plot(z1)
xlabel('t(seconds)')
ylabel('y(t)')
figure('Name','Case2 No Filter','NumberTitle','off')
plot(z2)
xlabel('t(seconds)')
ylabel('y(t)')
ylim([-2 2])
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/100)
figure('Name','Case3 Triangle','NumberTitle','off')
plot(z3)
xlabel('t(seconds)')
ylabel('y(t)')
xt = get(gca, 'XTick');                                 % 'XTick' Values
set(gca, 'XTick', xt, 'XTickLabel', xt/100)
%-----------------------------------------------------
##% plot output of the recieve filters --with different way 
% ------impulse reponses ------
h1 = g; % this is a matched filter with unit energy 

%----assume no noise -------------------------------
x = SS;
%----first case --------
out1 = conv(h1, x);
z1=conv(x,h1); %applying the recieved signal to matched filter
z1=sign(z1(T:T:end)); %sampling at T & using thresholding operation
recieved_bits1=(z1+1)/2; %recovering the sent bits
[number1,ratio1] = biterr(sent_bits,recieved_bits1)
##%----second case -------
%-----no conv here --------
z2= x; 
z2=sign(z2(T:T:end)); %sampling at T & using thresholding operation
recieved_bits2=(z2+1)/2; %recovering the sent bits
[number2,ratio2] = biterr(sent_bits,recieved_bits2)

%------thrid case ----------------
z3 = z3(fs:fs:end);
z3 = sign(z3);
z3 = z3(1:end-1);
recieved_bits3=(z3+1)/2; %recovering the sent bits
[number3,ratio3] = biterr(sent_bits,recieved_bits3)
  
figure('Name','Case 1 ,Matched Filter output','NumberTitle','off')
plot(out1)
xlabel('t(seconds)')
ylabel('y(t)')
#ylim([-2 2])
figure('Name','Case 2 ,NO Filter ','NumberTitle','off')
stairs([SS,SS(end)],'r','linewidth',2);
xlabel('t(seconds)')
ylabel('y(t)')
%---------------------------------------------------------
%---------------------------------------------------------
% --------------plot Eb/No vs bit error rate -------------
%---------------------------------------------------------
%---------------------------------------------------------
##
EbNoVec = (-10:2:20)';      % Eb/No values (dB)
berEst1 = zeros(size(EbNoVec)); % for first case 
berEst2 = zeros(size(EbNoVec)); % for second case 
berEst3 = zeros(size(EbNoVec)); % for third case 
PE = zeros(size(EbNoVec));

%------for bit error rate --------
%----generate too many bits ------
N = 10^4; %number of bits to be sent 
%data generation
sent_bits = randi([0 1],1,N);
%---if sent_bit=0 make it -1 otherwise =1
S = 2*sent_bits-1;
fs = 100;
S = upsample(S,fs);
h1 = ones(1,fs)./fs;
tmp = conv(S,h1);
sqrFilter = ones(1,fs)./fs;
triFilter = linspace(0,1,fs)./fs;
noFilter = [1 zeros(1,fs-1)]./fs;


for n = 1:length(EbNoVec)
  Eb=1; %energy per symbol 
  Eb_No = 10.^(0.1*EbNoVec(n));
  No = Eb/Eb_No;
  var = No/2;
  noise = sqrt(var).*randn(1,length(tmp)); %noise 
  x = tmp+noise; % add awgn noise to the sent_bits 
  %----first case --------
  z1=conv(x,sqrFilter); %applying the recieved signal to matched filter
  z1 = downsample(z1,fs);
  z1 = z1(2:end-1);
  z1 = sign(z1);
  recieved_bits1=(z1+1)/2; %recovering the sent bits
  [number1,ratio1] = biterr(sent_bits,recieved_bits1);
  %--second case ----------
  z2=conv(x,noFilter); %applying the recieved signal to matched filter
  z2 = downsample(z2,fs);
  %z1=sign(z1(T:T:end)); %sampling at T & using thresholding operation = 0
  z2 = z2(1:end-2);
  z2 = sign(z2);
  recieved_bits2=(z2+1)/2; %recovering the sent bits
  [number2,ratio2] = biterr(sent_bits,recieved_bits2);
  %---third case -----------
  z3=conv(x,triFilter); %applying the recieved signal to matched filter
  z3 = z3(fs:fs:end);
  z3 = sign(z3);
  z3 = z3(1:end-1);
  recieved_bits3=(z3+1)/2; %recovering the sent bits
  [number3,ratio3] = biterr(sent_bits,recieved_bits3);
  
  berEst1(n) = ratio1;
  berEst2(n) = ratio2;
  berEst3(n) = ratio3;
  
  %-------probability of error ---------------
  %---probability of error will be same of all cases 
  %---since E = 1 in all three cases --------------
  PE(n) = (0.5)*erfc(sqrt(Eb_No));
endfor
EbNoVec
%-----------Probabilit of error in different way---------
SNR=10.^(EbNoVec/10);
P_theory=(0.5)*erfc(sqrt(2.*SNR)./sqrt(2));
%------------------------------------------
figure('Name','BER Vs E/No','NumberTitle','off')
semilogy(EbNoVec,berEst1,'-')
hold on 
semilogy(EbNoVec,berEst2,'-')
hold on
semilogy(EbNoVec,berEst3,'-')
grid
legend('BER case1 matched filter','BER case2 no filter','BER case3 triangle')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')