
% The parameters
R1 = 1;
R2 = 2;
c = 0.25;
R3 = 10; %% R3 is 10 ohms
L = 0.2;
R4 = 0.1;
RO = 1000;

% The matrices for C and G, it is easier to do each row and fill with zeros instead of writing out the zeros
G=zeros(6);
C=zeros(6);
G(1,:)=[1 0 0 0 0 0]; 
C(1,:)=[0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1]; 
C(2,:)=[-c +c 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1]; 
C(3,:)=[0 0 0 0 0 0]; 
G(4,:)=[0 0 -1*100/R3 1 0 0]; 
C(4,:)=[0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0]; 
C(5,:)=[0 0 0 0 0 0];
G(6,:)=[0 -1 1 0 0 0]; 
C(6,:)=[0 0 0 0 0 L]; 


matrix=zeros(3,20);
matrix(1,:)=linspace(-10,10,20);
for i=1:20
    F=[matrix(1,i);
        0;
        0;
        0;
        0;
        0];
    V=G\F;
    matrix(2,i)=V(5);
    matrix(3,i)=V(3);
end

figure(1)
subplot(2,1,1);
plot(matrix(1,:),matrix(2,:));
title('Voltage as a Function of Input Voltage VO');
ylabel('Voltage (V)');
xlabel('V input (V)');
subplot(2,1,2);
plot(matrix(1,:),matrix(3,:));
title('Voltage as a Function of Input Voltage V3');
ylabel('Voltage (V)');
xlabel('V input (V)');



matrix2=zeros(2,2000);
matrix2(1,:)=linspace(0,500,2000);
Vin=1;
for i=1:2000
    F=[Vin;
        0;
        0;
        0;
        0;
        0];
    V=(G+1j*matrix2(1,i)*C)\F;
    matrix2(2,i)=V(5);
end

figure(2)
plot(matrix2(1,:),real(matrix2(2,:)));
title('Output Voltage as a Function of Angular Frequency');
ylabel('Output voltage');
xlabel('omega');

figure(3)
semilogx(matrix2(1,:),20*log10(real(matrix2(2,:))./Vin));
title('Gain as a Function of Angular Frequency');
ylabel('Gain (dB)')
xlabel('omega');


matrix3=zeros(1,10000);
for i=1:10000
    C(2,:)=[-(c+randn()*0.05) c+randn()*0.05 0 0 0 0];
    F=[Vin;
        0;
        0;
        0;
        0;
        0];
    V=(G+1j*pi*C)\F;
    matrix3(1,i)=V(5);
end

figure(4);
hist(real(matrix3),50);
title('Gain for perturbations in C');
ylabel('Counter')
xlabel('Gain dB');

C(2,:)=[-c +c 0 0 0 0];
stepsize=1/1000;
freq=(-(1000)/2:(1000)/2)*(1/1/1000/1000+1);
matrix3=zeros(3,3,1000+1);
oldVoltage=[0; 0; 0; 0; 0; 0];
input_sim=zeros(3,1000+1);
for type=1:3
    for i=1:1000
        if(type==1)
            if(i*stepsize<0.03)
                input_sim(type,i)=0;
            else
                input_sim(type,i)=1;
            end
        elseif(type==2)
            input_sim(type,i)=sin(2*pi*1/0.03*i*stepsize);
        else
            input_sim(type,i)=exp(-(i*stepsize-0.1)^2/(2*0.03^2));
        end
    end
end

for type=1:3
    oldVoltage=[0; 0; 0; 0; 0; 0];
    for i=1:1000
        F=[input_sim(type,i); 0; 0; 0; 0; 0];
        matrix3(type,1,i+1)=i*stepsize;
        matrix3(type,2,i+1)=input_sim(type,i);
        V=(C/stepsize+G)\(C*oldVoltage/stepsize+F);
        matrix3(type,3,i+1)=V(5);
        oldVoltage=V;
    end
    
end

matrix4=zeros(3,1000+1);
matrix4(:,:)=matrix3(1,:,:);
figure(5)
hold on;
subplot(2,1,1);
plot(matrix4(1,:),matrix4(2,:));
title('Figure 5: Voltages for Step Input (V1)');
ylabel('Voltage (V)');
xlabel('time (s)');
subplot(2,1,2);
plot(matrix4(1,:),matrix4(3,:));
title('Figure 5: Voltages for Step Input (VO)');
ylabel('Voltage (V)');
xlabel('time (s)');
hold off;


matrix5=zeros(3,1000+1);
matrix5(:,:)=matrix3(2,:,:);
figure(6)
hold on;
subplot(2,1,1);
plot(matrix5(1,:),matrix5(2,:));
title('Voltages for Sinusoidal Input with (VI)');
ylabel('Voltage (V)');
xlabel('time (s)');
subplot(2,1,2);
plot(matrix5(1,:),matrix5(3,:));
title('Voltages for Sinusoidal Input with (VO)');
ylabel('Voltage (V)');
xlabel('time (s)');
hold off;



matrix6=zeros(3,1000+1);
matrix6(:,:)=matrix3(3,:,:);
figure(7)
hold on;
subplot(2,1,1);
plot(matrix6(1,:),matrix6(2,:));
title('Figure 7: Voltages for Gaussian Pulse Input (VI)');
ylabel('Voltage (V)');
xlabel('time (s)');
subplot(2,1,2);
plot(matrix6(1,:),matrix6(3,:));
title('Figure 7: Voltages for Gaussian Pulse Input (VO)');
ylabel('Voltage (V)');
xlabel('time (s)');
hold off;



figure(8)
X=fft(matrix4(2,:));
Y=fft(matrix4(3,:));
hold on;
plot(freq,fftshift(abs(X)));
plot(freq,fftshift(abs(Y)));
hold off;
title('Figure 8: Frequency Domain for Step Input');
ylabel('Magnitude');
xlabel('frequency (Hz)');

figure(9)
X=fft(matrix5(2,:));
Y=fft(matrix5(3,:));
hold on;
plot(freq,fftshift(abs(X)));
plot(freq,fftshift(abs(Y)));
hold off;
title('Figure 9: Frequency Domain for Sinusoidal Input with f=33.3Hz');
ylabel('Magnitude');
xlabel('frequency (Hz)');

figure(10)
X=fft(matrix6(2,:));
Y=fft(matrix6(3,:));
hold on;
plot(freq,fftshift(abs(X)));
plot(freq,fftshift(abs(Y)));
hold off;
title('Figure 10: Frequency Domain for Gaussian Pulse Input');
ylabel('Magnitude');
xlabel('frequency (Hz)');
