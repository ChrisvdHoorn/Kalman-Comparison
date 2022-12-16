    clear,clc

%% system
dt = 0.1;
time = 0:dt:100;

A = [1 dt;
     0 1 ];
B = [0; dt];
C = [0 1];

Q = diag([0.01 0.01]);
R = 0.01;

%% data
dat_len = length(time);

u = linspace(0,0.1,dat_len) + sqrt(Q(1,1))*randn(dat_len,1)';
x = zeros(2,dat_len);

for i = 2:dat_len
    noise = [sqrt(Q(1,1))*randn(1,1);
             sqrt(Q(2,2))*randn(1,1)];
    x(:,i) = A*x(:,i-1) + B*u(i-1) + noise;
end

%% static kalman gain
if rank(obsv(A,C)) < size(A,1)
    disp("system is unobservable");
    Kstat = [0;0];
else
    [~,Kstat,~,~] = idare(A',C',Q,R,[],[]);
    Kstat = Kstat';
end

%% filter and dynamic gain
static = zeros(2,dat_len);
dynamic = zeros(2,dat_len);

P = zeros(2,2);

for i = 2:dat_len
    meas_noise = sqrt(R)*randn(1,1);
% static
    static(:,i) = A*static(:,i-1) + B*u(i-1) + Kstat*(C*x(:,i) + meas_noise  - C*static(:,i-1));
% dynamic
    % Predict
    x_int_est = A*dynamic(:,i-1) + B*u(i-1);
    P = A*P*A' + Q;

    % Update
    y_est = (C*x(:,i) + meas_noise) - C*x_int_est;
    S = C*P*C' + R;
    Kdyn = P*C'*inv(S); %#ok<MINV>   
    dynamic(:,i) = x_int_est + Kdyn*y_est;
    P = (eye(2) - Kdyn*C)*P;
end


%% plot

figure(2),clf
subplot(2,1,1),hold on,grid on;
    ylabel('x_1')
    title(['Filtered states using C = [' num2str(C(1)) ', ' num2str(C(2)) ']'])
    plot(time,x(1,:),'DisplayName','data')
    plot(time,static(1,:),'--','DisplayName','static')
    plot(time,dynamic(1,:),'-.','DisplayName','dynamic')
    legend();
subplot(2,1,2),hold on,grid on;
    ylabel('x_2')
    plot(time,x(2,:),'DisplayName','data')
    plot(time,static(2,:),'--','DisplayName','static')
    plot(time,dynamic(2,:),'-.','DisplayName','dynamic')
    xlabel('time(s)')
    legend();