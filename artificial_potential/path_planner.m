%%%%%
% Main entrypoint for FINCH ADCS artificial potential path planner dev tool
% Author: Alexander Wlodarczyk
%%%%%
clear all;

%% Spacecraft specifications
exc_str_h = deg2rad(40); % star tacker exclusion cone half-angle
exc_cam_h = deg2rad(15); % camera exclusion cone half-angle
u_sa = [0; 1; 0]; % primary solar array vector in body frame
u_str = [0; -1; 0]; % star tracker vector in body frame
u_cam = [0; 0; 1]; % payload camera vector in body frame
refs = [u_cam'; u_str'; u_str']; % list of reference vectors associated with constraints

%% Scene Setup
% Sun and moon vectors
v_sun = rand(1,3); % random vector to Sun in ECI
u_sun = v_sun / norm(v_sun); % unit vector to Sun in ECI
v_moon = rand(1,3); % random vector to Moon in ECI
u_moon = v_moon / norm(v_moon); % unit vector to moon in ECI
obs = [u_sun; u_sun; u_moon]'; % list of constraints

% Initial condition
q0 = random_quat(); % Random quaternion
xi0 = q2xi(q0); % Lie algebra space
v0 = quat2rotm(q0'); % Direction Cosine Matrix
% Ensure IC is not within exclusion zone
while dist_angle(v0(1,:), u_sun) <= exc_cam_h | dist_angle(v0(2,:), u_sun) <= exc_str_h
    q0 = random_quat();
    xi0 = q2xi(q0);
    v0 = quat2rotm(q0');
end

% Final condition
qf = random_quat(); % Random quaternion
xif = q2xi(qf); % Lie algebra space
vf = quat2rotm(qf'); % Direction Cosine Matrix
% Ensure FC is not wihtin exclusion zone
while dist_angle(vf(1,:), u_sun) <= exc_cam_h | dist_angle(vf(2,:), u_sun) <= exc_str_h
    qf = random_quat();
    xif = q2xi(qf);
    vf = quat2rotm(qf');
end


%% Run iteratively
% Configuration
step = 0.001;
max_iters = 1000;
max_step = 0.08;
tol = 0.01;
% Initialize trajectory
q = q0;
Q_traj = zeros(4, max_iters);
warnings = zeros(3, max_iters);
R_traj = zeros(3,3,max_iters);
E = zeros(1, max_iters);
U = zeros(1, max_iters);
for k = 1:max_iters
    Q_traj(:,k) = q;
    R = quat2rotm(q'); R_traj(:,:,k) = R;

    %[U(k), g] = nav(q, qf, refs, obs);
    U(k) = potential(q, qf, refs, obs);
    g = gradU(q, qf, refs, obs);

    u = -step*g;
    if norm(u) > max_step
        u = u / norm(u) * max_step;
    end

    dq = expm_so3(u);
    q = qnormalize(qmult(dq', q'));

    % orientation error
    E(k) = norm(qlogmap(q, qf));

    % restriction error
    warnings(1,k) = dist_angle(R(1,:), u_sun);
    warnings(2,k) = dist_angle(R(2,:), u_sun);
    warnings(3,k) = dist_angle(R(2,:), u_moon);
end

%% Visualize
% Angular distance between reference vector and target
figure;
plot(E);
xlabel('time step');
ylabel('orientation error');

% Euler angles
figure;
eul = quat2eul(Q_traj');
plot(eul);
legend('yaw','pitch','roll');
xlabel('time step');
ylabel('angle [rad]');

% Distances between sensitive instruments and forbidden regions
figure;
subplot(3,1,1);
hold on;
plot(warnings(1,:));
yline(exc_cam_h, 'r--');
hold off;
ylim([0, 2*pi]);

subplot(3,1,2);
hold on;
plot(warnings(2,:));
yline(exc_str_h, 'r--');
hold off;
ylim([0, 2*pi]);

subplot(3,1,3);
hold on;
plot(warnings(3,:));
yline(exc_str_h, '-.');
hold off;
ylim([0, 2*pi]);

% Total potential
figure;
plot(U);

% Animated orientation in R3
figure;
hold on;
dt = 0.01;
% static vectors
plot3([0, u_sun(1)], [0, u_sun(2)], [0, u_sun(3)], 'Color', '#ffa500','LineWidth',2, 'DisplayName', '$\vec{u}_{sun}$');
plot3([0, u_moon(1)], [0, u_moon(2)], [0, u_moon(3)], 'Color', '#808080','LineWidth',2, 'DisplayName', '$\vec{u}_{moon}$');
plot3([0, vf(1,1)], [0, vf(1,2)], [0, vf(1,3)], 'Color', 'k','LineWidth',2, 'DisplayName', '$\vec{q}_t$');
% initialize animation
hx = animatedline('Color','r','LineWidth',2,'DisplayName','$\vec{x}_B$');
hy = animatedline('Color','g','LineWidth',2,'DisplayName','$\vec{y}_B$');
hz = animatedline('Color','b','LineWidth',2,'DisplayName','$\vec{z}_B$');
hist_x = animatedline('Color', [1, 0, 0, 0.5], 'DisplayName', 'Traced path X');
hist_y = animatedline('Color', [0, 1, 0, 0.5], 'DisplayName', 'Traced path Y');
view(3); axis equal tight;
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]); zlim([-1.5, 1.5]);
xlabel('x'); ylabel('y'); zlabel('z');
legend('Interpreter', 'latex');
for k = 1:max_iters
    % draw spacecraft axes in frame R
    clearpoints(hx);
    clearpoints(hy);
    clearpoints(hz);
    addpoints(hx, 0, 0, 0); addpoints(hx, R_traj(1,1,k), R_traj(1,2,k), R_traj(1,3,k));
    addpoints(hy, 0, 0, 0); addpoints(hy, R_traj(2,1,k), R_traj(2,2,k), R_traj(2,3,k));
    addpoints(hz, 0, 0, 0); addpoints(hz, R_traj(3,1,k), R_traj(3,2,k), R_traj(3,3,k));
    addpoints(hist_x, R_traj(1,1,k), R_traj(1,2,k), R_traj(1,3,k));
    addpoints(hist_y, R_traj(2,1,k), R_traj(2,2,k), R_traj(2,3,k));

    drawnow;
    pause(dt);
end
hold off;