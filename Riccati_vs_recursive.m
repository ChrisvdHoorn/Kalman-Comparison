clear,clc
%% system
dt = 0.1;
A = [1     dt
     0     1 ];
B = [0;dt];
C = [1 0];
D = 0;

Q = diag([0.01 0.01]);
R = 0.01;

N = zeros(2,1);
G = zeros(2,1);
H = 0;
sys = ss(A, [B, G], C, [D, H], dt, 'inputname', 'u', 'outputname', 'y');
[kalmf,kalmG,~] = kalman(sys,Q,R,N);

%% static kalman gain
if rank(obsv(A,C)) < size(A,1)
    disp("system is unobservable");
    Kstat = [0;0];
else
    [~,Kstat,~,~] = idare(A',C',Q,R,[],[]);
    Kstat = Kstat';
end

%% find P_inf
iterations = 60;

Kdyn_array = zeros(2,iterations);
P = zeros(2,2);

for i = 2:iterations
    P = A*P*A' + Q;
    S = C*P*C' + R;
    Kdyn = P*C'*inv(S); %#ok<MINV>  
    P = (eye(2) - Kdyn*C)*P;

    Kdyn_array(:,i) = Kdyn;
end

%% plot
figure(1),clf;
subplot(2,1,1),hold on,grid on;
    plot([1,iterations],[Kstat(1) Kstat(1)])
    plot(Kdyn_array(1,:),'--')
    legend('K-static','K-dynamic')
    title('Static vs Recursive Kalman gain')
    ylabel('K_1')
subplot(2,1,2),hold on,grid on;
    plot([1,iterations],[Kstat(2) Kstat(2)])
    plot(Kdyn_array(2,:),'--')
    legend('K-static','K-dynamic')
    ylabel('K_2')
    xlabel('iteration')
