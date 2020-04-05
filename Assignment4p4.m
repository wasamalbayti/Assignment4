
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
RO = 1000;
C1 = 0.25;
L1 = 0.2;
alpha = 100;
% new values for beta and gamma
beta=10e6;
gamma=1e9;

inputVoltage=1;
G=zeros(6);
C=zeros(6);
G(1,:)=[1 0 0 0 0 0]; % V1
C(1,:)=[0 0 0 0 0 0]; % V1
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1]; 
C(2,:)=[-C1 +C1 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1]; 
C(3,:)=[0 0 0 0 0 0]; 
G(4,:)=[0 0 -1*alpha/R3 1 0 0]; 
C(4,:)=[0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0]; 
C(5,:)=[0 0 0 0 0 0];
G(6,:)=[0 -1 1 0 0 0]; 
C(6,:)=[0 0 0 0 0 L1]; 

stepN=1000;
timestep=1/stepN;

A=(C./timestep+G);

f=1/0.03;

V(1:6,1)=[0;0;0;0;0;0];
oldVolt=V;

for step=1:stepN
    t=step*timestep;
    inputVoltage=sin(2*pi*f*t);
    F=[inputVoltage; 0; 0; 0; 0; 0];
    counter=0;
    while(counter<150)
        I3=V(6);
        voltageD=-2*beta*I3-3*gamma*I3.^2;
    
        J=[0,0,0,0,0,0;
           0,0,0,0,0,0;
           0,0,0,0,0,0;
           0,0,voltageD,0,0,0;
           0,0,0,0,0,0;
           0,0,0,0,0,0];
       
        B=[0;
            0;
            0;
            beta*I3.^2+gamma*I3.^3;
            0;
            0;];
        voltageFreq=(C./timestep+G)*V-(C./timestep)*oldVolt-F-B;
        H=C./timestep+G-J;
        dV=-inv(H)*voltageFreq;
        V=V+dV;
        if(max(abs(dV))<5e-4)   
            break;     
        end 
        counter=counter+1;
    end
    inputvoltage_1(step)=V(1,1);
    outputVoltage_1(step)=V(5,1);
    oldVolt=V;
    
end
figure(3)
plot(linspace(0,1,stepN),inputvoltage_1)
title('Input Voltage: Step Response')
xlabel('Time (s)')
ylabel('Voltage (V)')

figure(4)
plot(linspace(0,1,stepN),outputVoltage_1)
title('Output Voltage: Step Response')
xlabel('Time (s)')
ylabel('Voltage (V)')

inputfft = fft(inputvoltage_1);
outputfft = fft(outputVoltage_1);
figure(5)
plot(linspace(-1/timestep*0.5,1/timestep*0.5,length(inputfft)),fftshift(abs(inputfft)))
title('Input Frequency Content: Step Response')
xlabel('Frequency (Hz)')
ylabel('Power')

figure(6)
plot(linspace(-1/timestep*0.5,1/timestep*0.5,length(outputfft)),fftshift(abs(outputfft)))
title('Output Frequency Content Step Response')
xlabel('Frequency (Hz)')
ylabel('Power')



