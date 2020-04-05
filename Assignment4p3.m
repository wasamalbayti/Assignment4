% paramaters
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
RO = 1000;
c = 0.25;
L = 0.2;
Cn = 0.00001;

%matrices
Vin=1;
G=zeros(6);
C=zeros(6);
G(1,:)=[1 0 0 0 0 0]; 
C(1,:)=[0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1]; 
C(2,:)=[-c +c 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1]; 
C(3,:)=[0 0 Cn 0 0 0]; 
G(4,:)=[0 0 -1*100/R3 1 0 0]; 
C(4,:)=[0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0]; 
C(5,:)=[0 0 0 0 0 0];
G(6,:)=[0 -1 1 0 0 0]; 
C(6,:)=[0 0 0 0 0 L]; 



C(2,:)=[-0.25 +0.25 0 0 0 0];
timesteps=1000;
Step=1/1000;
fD_Q2=(-(timesteps)/2:(timesteps)/2)*(1/Step/timesteps+1);
simInput=zeros(3,timesteps+1);
for x=1:3
    for i=1:timesteps
        simInput(x,i)=exp(-(i*Step-0.1)^2/(2*0.03^2));
    end
end
simdata=zeros(3,3,timesteps+1);

for x=1:3

    if(x==1)
        cn1=Cn;
    elseif(x==2)
        cn1=Cn*100;
    else
        cn1=Cn*1000;
    end
    C(3,:)=[0 0 cn1 0 0 0]; 
    oldVoltage=[0;
        0;
        0;
        0;
        0;
        0];
    for i=1:timesteps  
        Vin=simInput(x,i);
        F=[Vin;
            0; 
            -0.001*randn();
            0;
            0;
            0];
        simdata(x,1,i+1)=i*Step;
        simdata(x,2,i+1)=Vin;
        V=(C/Step+G)\(C*oldVoltage/Step+F);
        simdata(x,3,i+1)=V(5);
        oldVoltage=V;
    end
end
C(3,:)=[0 0 Cn 0 0 0]; 

simdata_1=zeros(3,timesteps+1);
simdata_1(:,:)=simdata(1,:,:);
figure(1)
hold on;
plot(simdata_1(1,:),simdata_1(2,:));
plot(simdata_1(1,:),simdata_1(3,:));
hold off;
title('Voltages for Gaussian Pulse Input with Noise');
ylabel('Voltage (V)');
xlabel('time (s)');

figure(2)
X=fft(simdata_1(2,:));
Y=fft(simdata_1(3,:));
hold on;
plot(fD_Q2,fftshift(abs(X)));
plot(fD_Q2,fftshift(abs(Y)));
hold off;
title('Frequency Domain for Gaussian Pulse Input with Noise');
ylabel('Magnitude');
xlabel('frequency (Hz)');
simdata_1=zeros(3,timesteps+1);
simdata_1(:,:)=simdata(1,:,:);
figure(3)
subplot(3,1,1)
hold on;
plot(simdata_1(1,:),simdata_1(2,:));
plot(simdata_1(1,:),simdata_1(3,:));
hold off;
title({'Voltages with Noise for Different Cn'});
ylabel('Voltage (V)');
xlabel('time (s)');

simdata_2=zeros(3,timesteps+1);
simdata_2(:,:)=simdata(2,:,:);
subplot(3,1,2)
hold on;
plot(simdata_2(1,:),simdata_2(2,:));
plot(simdata_2(1,:),simdata_2(3,:));
hold off;
title('Cn = 1mF');
ylabel('Voltage (V)');
xlabel('time (s)');

simdata_3=zeros(3,timesteps+1);
simdata_3(:,:)=simdata(3,:,:);
subplot(3,1,3)
hold on;
plot(simdata_3(1,:),simdata_3(2,:));
plot(simdata_3(1,:),simdata_3(3,:));
hold off;
title('Cn=10mF');
ylabel('Voltage (V)');
xlabel('time (s)');


C(2,:)=[-c +c 0 0 0 0];
timesteps_ALT=10000;
nStep=1/timesteps_ALT;
numT4=timesteps_ALT+1;
fD_Q3=(-(numT4-1)/2:(numT4-1)/2)*(1/nStep/numT4);
datasimulation1=zeros(3,3,timesteps_ALT+1);
oldVoltage=[1; 0; 0; 0; 0; 0];
simInput_ALT=zeros(3,timesteps_ALT+1);
freq=1/0.03;
for x=1:1
    for i=1:timesteps_ALT
        simInput_ALT(x,i)=exp(-(i*nStep-0.1)^2/(2*0.03^2));
    end
end

for x=1:1
    C(3,:)=[0 0 Cn 0 0 0]; 

    oldVoltage=[0; 0; 0; 0; 0; 0];

    for i=1:timesteps_ALT
        Vin=simInput_ALT(x,i);
        In=0.001*randn();
        F=[Vin; 0; -In; 0; 0; 0];
        datasimulation1(x,1,i+1)=i*nStep;
        datasimulation1(x,2,i+1)=Vin;
        A=C/nStep+G;
        V=(A)\(C*oldVoltage/nStep+F);
        datasimulation1(x,3,i+1)=V(5);
        oldVoltage=V;
    end
    
end
C(3,:)=[0 0 Cn 0 0 0]; 
figure(4)
subplot(2,1,1)
simdata_1=zeros(3,timesteps+1);
simdata_1(:,:)=simdata(1,:,:);
hold on;
plot(simdata_1(1,:),simdata_1(2,:));
plot(simdata_1(1,:),simdata_1(3,:));
hold off;
title({'Voltages with Noise for Different Timesteps'});
ylabel('Voltage (V)');
xlabel('time (s)');

subplot(2,1,2)
simdata_1_ALT=zeros(3,timesteps_ALT+1);
simdata_1_ALT(:,:)=datasimulation1(1,:,:);
hold on;
plot(simdata_1_ALT(1,:),simdata_1_ALT(2,:));
plot(simdata_1_ALT(1,:),simdata_1_ALT(3,:));
hold off;
title('10000 timesteps');
ylabel('V (V)');
xlabel('time (s)');

