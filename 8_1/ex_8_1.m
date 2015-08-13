%% �^�l�̍쐬
N = 200;
z = zeros(N,1);
d = randn(N,1);
z(1) = 10;
for k=2:N
	z(k) = z(k-1) + 2*cos(0.05*k) + d(k-1);
end

%% �ϑ��l�̐���
% �Z���T�덷�̐ݒ�
a = 0.75;
A = [a 0 0; 0 1 0; 0 0 0];
B = [1-a 0; 0 0; 0 1];
Q = diag([60,60]);
R = 1e-4;
C = [1 1 -1]';
% �Z���T�덷�̐ݒ�
v = randn(N,2) * sqrtm(Q);
x = zeros(N,3);
x(1,:) = [0;20;0];
for k=2:N
	x(k,:) = A*x(k-1,:)' + B*v(k-1,:)';
end

e1 = x(:,1) + x(:,2);
e2 = x(:,3);

% �Z���T�o��
y1 = z + e1;
y2 = z + e2;

% ����t�B���^���p����o��
y = y1 - y2;

%% ����덷e1,e2�̐U���X�y�N�g���̌v�Z
fs = 1;
m = length(e1);
n = pow2(nextpow2(m));
f = (0:n-1) * (fs/n);
ye1 = fft(e1,n);
ye2 = fft(e2,n);
me1 = abs(ye1);
me2 = abs(ye2);

%% ����t�B���^�ɂ��덷����
xhat = zeros(N,3);

% ��ԏ����l
xhat(1,:) = x(1,:)' + [10;10;20];
P = 1000 * eye(3);

% ���ԍX�V
for k=2:N
	[xhat(k,:),P] = kf(A,B,0,C,Q,R,0,y(k),xhat(k-1,:),P);
end
error = xhat - x;

% ����l��p���Ċe�Z���T�o�͂�␳
yhat1 = y1 - xhat(:,1) - xhat(:,2);
yhat2 = y2 - xhat(:,3);

%% ���ʂ̐}��
figure(1), clf
plot(z), xlabel('k'), ylabel('z')

figure(2), clf
plot(e1), xlabel('k'), ylabel('e1')
figure(3), clf
plot(e2), xlabel('k'), ylabel('e2')

figure(4), clf
plot(f, me1), xlabel('Normalized freq'), ylabel('e1')
axis([0 0.5 0 1000])
figure(5), clf
plot(f, me2), xlabel('Normalized freq'), ylabel('e2')
axis([0 0.5 0 1000])

figure(6), clf
plot(1:N, y1, 'k:', 1:N, z, 'r--')
legend('Sensor1', 'True');
figure(7), clf
plot(1:N, y2, 'k:', 1:N, z, 'r--')
legend('Sensor2', 'True');

figure(8), clf
plot(1:N, z, 'k:', 1:N, y1, 'r--', 1:N, yhat1, 'b--')
xlabel('k'), ylabel('z')
title('Mesurement with sensor 1')
legend('True', 'Raw', 'Corrected')

figure(9), clf
plot(1:N, z, 'k:', 1:N, y2, 'r--', 1:N, yhat2, 'b--')
xlabel('k'), ylabel('z')
title('Mesurement with sensor 2')
legend('True', 'Raw', 'Corrected')
