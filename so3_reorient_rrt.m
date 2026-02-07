function result = so3_reorient_rrt()
% SO(3) satellite reorientation planner with a Sun keep-out constraint.
% - Attitude represented as unit quaternion q = [qw qx qy qz] (scalar first).
% - Sensitive body axis a_sens_body (unit 3x1) must avoid Sun direction s_inertial (unit 3x1)
%   by at least alpha (radians).
% - Uses RRT in quaternion space with geodesic metric and SLERP edge checking.

%% -------------------- User inputs (example) --------------------
% Start and goal attitudes (unit quaternions, scalar first)
% Sun direction in inertial frame (unit vector)
sInertial = [1; 0; 0]; sInertial = sInertial / norm(sInertial);

% Sensitive body axis (unit vector, in body frame)
a1 = [0;0;1];   % payload boresight
a2 = [1;0;0];   % e.g., star tracker axis
ASensBody = [a1/norm(a1), a2/norm(a2)];   % 3x2

alpha = deg2rad([25, 40]);  % 1x2


% Initial condition
qStart = randUnitQuat(); % Random quaternion
if ~isQuatFeasible(qStart, sInertial, ASensBody, alpha)
    qStart = randUnitQuat();
end

% Final condition
qGoal = randUnitQuat(); % Random quaternion
if ~isQuatFeasible(qGoal, sInertial, ASensBody, alpha)
    qGoal = randUnitQuat();
end

% qStart = quatNormalize(axang2quat([0 1 0 pi/3]));              % identity
% qGoal  = quatNormalize(axang2quat([0 1 0 pi/3+pi/2])); % rotate +90deg about y

qStart = quatNormalize(qStart);
qGoal = quatNormalize(qGoal);

% RRT parameters
params.maxIters       = 4000;
params.goalBias       = 0.15;     % probability of sampling goal
params.stepAngle      = deg2rad(8);% max geodesic step per RRT extension
params.edgeCheckStep  = deg2rad(2);% discretization along each edge for collision checking
params.goalTol        = deg2rad(3);% success if within this geodesic angle to goal

rng(1); % repeatable

%% -------------------- Plan --------------------
[result.success, pathQ, info] = rrtSO3(qStart, qGoal, sInertial, ASensBody, alpha, params);

result.info = info;
result.pathQ = pathQ;

if result.success
    fprintf("SUCCESS. Nodes=%d, Iters=%d, PathLen=%d\n", info.numNodes, info.iters, size(pathQ,1));
    % Optional: visualize angle-to-sun constraint along path
    % plotConstraintAlongPath(pathQ, sInertial, ASensBody, alpha);
    animatePath(pathQ, sInertial, ASensBody, alpha, ...
    "fps", 30, "speed", 0.5, "loops", 2, "pingpong", true);
else
    fprintf("FAILED. Nodes=%d, Iters=%d\n", info.numNodes, info.iters);
end

end

%% ========================================================================
function [success, pathQ, info] = rrtSO3(qStart, qGoal, sInertial, aSensBody, alpha, params)
disp(qStart)
% Validate start/goal are feasible
if ~isQuatFeasible(qStart, sInertial, aSensBody, alpha)
    error("Start attitude violates Sun keep-out constraint.");
end
if ~isQuatFeasible(qGoal, sInertial, aSensBody, alpha)
    error("Goal attitude violates Sun keep-out constraint.");
end

% Tree storage
Q   = zeros(params.maxIters+1, 4);   % node quaternions
par = zeros(params.maxIters+1, 1);   % parent indices
Q(1,:) = qStart;
par(1) = 0;
nNodes = 1;

success = false;
goalIdx = -1;

for it = 1:params.maxIters
    % Sample
    if rand < params.goalBias
        qRand = qGoal;
    else
        qRand = randUnitQuat();
    end

    % Nearest neighbor in geodesic angle
    idxNear = nearestQuat(Q(1:nNodes,:), qRand);
    qNear = Q(idxNear,:);

    % Steer: take a bounded geodesic step toward qRand
    qNew = steerQuat(qNear, qRand, params.stepAngle);

    % Feasibility checks: node + edge
    if ~isQuatFeasible(qNew, sInertial, aSensBody, alpha)
        continue;
    end
    if ~edgeFeasible(qNear, qNew, sInertial, aSensBody, alpha, params.edgeCheckStep)
        continue;
    end

    % Add node
    nNodes = nNodes + 1;
    Q(nNodes,:) = qNew;
    par(nNodes) = idxNear;

    % Check goal
    if quatGeoAngle(qNew, qGoal) <= params.goalTol
        % Also require connecting edge to goal feasible
        if edgeFeasible(qNew, qGoal, sInertial, aSensBody, alpha, params.edgeCheckStep)
            nNodes = nNodes + 1;
            Q(nNodes,:) = qGoal;
            par(nNodes) = nNodes-1;
            goalIdx = nNodes;
            success = true;
            info.iters = it;
            break;
        end
    end
end

info.numNodes = nNodes;
if ~isfield(info,'iters'); info.iters = params.maxIters; end

if ~success
    pathQ = zeros(0,4);
    return;
end

% Reconstruct path
idx = goalIdx;
pathIdx = idx;
while par(idx) ~= 0
    idx = par(idx);
    pathIdx(end+1) = idx; %#ok<AGROW>
end
pathIdx = fliplr(pathIdx);

pathQ = Q(pathIdx,:);

% Optional: shortcut smoothing (simple)
pathQ = shortcutPath(pathQ, sInertial, aSensBody, alpha, params.edgeCheckStep, 120);

end

%% ========================================================================
% ------------------------ Core geometry utilities ------------------------

function q = quatNormalize(q)
q = q(:).';
q = q / norm(q);
% Canonicalize sign (optional): keep qw >= 0 for consistency
if q(1) < 0, q = -q; end
end

function q = randUnitQuat()
% Uniform random unit quaternion (Shoemake method)
u1 = rand; u2 = rand; u3 = rand;
q = [ ...
    sqrt(1-u1)*sin(2*pi*u2), ...
    sqrt(1-u1)*cos(2*pi*u2), ...
    sqrt(u1)*sin(2*pi*u3), ...
    sqrt(u1)*cos(2*pi*u3)];
% This gives [qx qy qz qw]; convert to [qw qx qy qz]
q = [q(4) q(1) q(2) q(3)];
q = quatNormalize(q);
end

function ang = quatGeoAngle(q1, q2)
% Geodesic angle on SO(3): theta in [0, pi]
q1 = quatNormalize(q1); q2 = quatNormalize(q2);
d = abs(dot(q1,q2));           % account for double cover
d = min(1, max(-1, d));
ang = 2 * acos(d);
end

function idx = nearestQuat(Q, qTarget)
% Find nearest neighbor in geodesic angle
dots = abs(Q*qTarget(:));      % Nx1
dots = min(1, max(-1, dots));
angs = 2*acos(dots);
[~, idx] = min(angs);
end

function qNew = steerQuat(qFrom, qToward, maxStep)
% Move from qFrom toward qToward by at most maxStep along the geodesic (SLERP)
theta = quatGeoAngle(qFrom, qToward);
if theta < 1e-12
    qNew = quatNormalize(qToward);
    return;
end
t = min(1, maxStep/theta);
qNew = slerp(qFrom, qToward, t);
end

function q = slerp(q1, q2, t)
% Spherical linear interpolation between unit quaternions
q1 = quatNormalize(q1); q2 = quatNormalize(q2);

% Shortest path handling
if dot(q1,q2) < 0
    q2 = -q2;
end
dot12 = dot(q1,q2);
dot12 = min(1, max(-1, dot12));

if dot12 > 0.9995
    % Nearly identical: linear interpolate then renormalize
    q = quatNormalize((1-t)*q1 + t*q2);
    return;
end

omega = acos(dot12);
sinom = sin(omega);
q = (sin((1-t)*omega)/sinom)*q1 + (sin(t*omega)/sinom)*q2;
q = quatNormalize(q);
end

function R = quat2rotm_scalarFirst(q)
% Convert quaternion [qw qx qy qz] to rotation matrix
q = quatNormalize(q);
qw=q(1); qx=q(2); qy=q(3); qz=q(4);

R = [ ...
    1-2*(qy^2+qz^2),   2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
    2*(qx*qy+qz*qw),   1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
    2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)];
end

%% ========================================================================
% ------------------------ Constraint checking ----------------------------

function ok = isQuatFeasible(q, sInertial, ASensBody, alpha)
% ASensBody: 3xK body axes (columns), each unit
% alpha: scalar or 1xK keep-out half-angles (radians)

R = quat2rotm_scalarFirst(q);           % body->inertial
Ain = R * ASensBody;                    % 3xK axes in inertial

s = sInertial(:) / norm(sInertial);
c = sum(Ain .* s, 1);                   % 1xK dot products
c = min(1, max(-1, c));
ang = acos(c);                          % 1xK angles

if isscalar(alpha)
    ok = all(ang >= alpha);
else
    alpha = alpha(:).';                 % force 1xK
    ok = all(ang >= alpha);
end
end


function ok = edgeFeasible(q1, q2, sInertial, ASensBody, alpha, stepAng)
% Sample along SLERP and ensure all intermediate attitudes feasible
theta = quatGeoAngle(q1, q2);
if theta < 1e-12
    ok = isQuatFeasible(q1, sInertial, ASensBody, alpha);
    return;
end
n = max(2, ceil(theta/stepAng)+1);
for k = 0:(n-1)
    t = k/(n-1);
    qk = slerp(q1, q2, t);
    if ~isQuatFeasible(qk, sInertial, ASensBody, alpha)
        ok = false;
        return;
    end
end
ok = true;
end

%% ========================================================================
% ------------------------ Path post-processing ---------------------------

function pathQ = shortcutPath(pathQ, sInertial, aSensBody, alpha, edgeCheckStep, nTrials)
% Simple random shortcutting: try to connect two random waypoints directly.
if size(pathQ,1) < 3
    return;
end
m = size(pathQ,1);
for t = 1:nTrials
    i = randi([1, m-1]);
    j = randi([i+1, m]);
    if j <= i+1, continue; end
    if edgeFeasible(pathQ(i,:), pathQ(j,:), sInertial, aSensBody, alpha, edgeCheckStep)
        % Replace segment i..j with direct connection
        pathQ = [pathQ(1:i,:); pathQ(j:end,:)];
        m = size(pathQ,1);
        if m < 3, break; end
    end
end
end

%% ========================================================================
% ------------------------ Visualization helper ---------------------------

function plotConstraintAlongPath(pathQ, sInertial, aSensBody, alpha)
angles = zeros(size(pathQ,1),1);
for i = 1:size(pathQ,1)
    R = quat2rotm_scalarFirst(pathQ(i,:));
    aIn = R*aSensBody;
    c = dot(aIn, sInertial); c = min(1,max(-1,c));
    angles(i) = acos(c);
end

figure; plot(rad2deg(angles), 'LineWidth', 1.5);
yline(rad2deg(alpha), '--', 'LineWidth', 1.2);
xlabel('Waypoint index'); ylabel('Angle(sensitive axis, Sun) [deg]');
title('Keep-out constraint along planned path');
grid on;
end

%% ========================================================================
% ------------------------ Axis-angle to quaternion -----------------------

function q = axang2quat(axang)
% axang = [ax ay az angle], axis need not be unit
ax = axang(1:3); th = axang(4);
ax = ax(:); ax = ax / norm(ax);
qw = cos(th/2);
qv = ax*sin(th/2);
q = quatNormalize([qw qv(:).']);
end

function animatePath(pathQ, sInertial, ASensBody, alpha, varargin)
% animatePath  Animate planned SO(3) attitude path with Sun keep-out cones.
%
% Inputs:
%   pathQ      : Mx4 unit quaternions [qw qx qy qz] (scalar-first), waypoint path
%   sInertial  : 3x1 Sun direction in inertial frame (unit or non-unit)
%   ASensBody  : 3xK (or Kx3) sensitive body axes (each column is an axis)
%   alpha      : keep-out half-angle in radians (scalar or 1xK)
%
% Options (name/value):
%   "fps"              (default 30)
%   "speed"            (default 1.0)  % 1=normal, 0.5=half speed, 2=double
%   "samplesPerSegment"(default 25)   % interpolation frames per waypoint segment
%   "showTriad"        (default true)
%   "loops"            (default 1)
%   "pingpong"         (default false) % if true, plays forward then backward
%   "saveGif"          (default false)
%   "gifFile"          (default "attitude_path.gif")
%   "sphereResolution" (default 30)
%
% Dependencies:
%   - quat2rotm_scalarFirst(q)  must return a 3x3 rotation matrix
%   - slerp(q1,q2,t)            quaternion slerp (unit quats, scalar-first)
%
% Convention assumed (important):
%   quat2rotm_scalarFirst(q) returns R such that v_inertial = R * v_body.

% -------------------- Parse options --------------------
p = inputParser;
addParameter(p, "fps", 30);
addParameter(p, "speed", 1.0);
addParameter(p, "samplesPerSegment", 25);
addParameter(p, "showTriad", true);
addParameter(p, "loops", 1);
addParameter(p, "pingpong", false);
addParameter(p, "saveGif", false);
addParameter(p, "gifFile", "attitude_path.gif");
addParameter(p, "sphereResolution", 30);
parse(p, varargin{:});
opt = p.Results;

if isempty(pathQ) || size(pathQ,2) ~= 4
    error("pathQ must be Mx4 quaternions.");
end

% -------------------- Normalize/shape inputs --------------------
sInertial = sInertial(:);
if numel(sInertial) ~= 3, error("sInertial must be 3x1."); end
sInertial = sInertial / norm(sInertial);

ASensBody = ap_ensure3xK(ASensBody);          % 3xK
% normalize each axis column
ASensBody = ASensBody ./ vecnorm(ASensBody);

K = size(ASensBody,2);

% alpha: scalar or 1xK
if isscalar(alpha)
    alphaVec = alpha * ones(1,K);
else
    alphaVec = alpha(:).';
    if numel(alphaVec) ~= K
        error("alpha must be scalar or have one value per sensitive axis (1xK).");
    end
end

% -------------------- Build dense path (SLERP samples) --------------------
[qDense, ~] = ap_densifyQuatPath(pathQ, opt.samplesPerSegment);
N = size(qDense,1);

% Precompute inertial vectors for animation
aSensIn = zeros(N,3,K);     % sensitive axes in inertial
xIn = zeros(N,3); yIn = zeros(N,3); zIn = zeros(N,3);

for i = 1:N
    R = quat2rotm_scalarFirst(qDense(i,:));   % 3x3, body->inertial
    Ain = R * ASensBody;                     % 3xK
    for j = 1:K
        aSensIn(i,:,j) = Ain(:,j).';
    end

    if opt.showTriad
        xIn(i,:) = (R*[1;0;0]).';
        yIn(i,:) = (R*[0;1;0]).';
        zIn(i,:) = (R*[0;0;1]).';
    end
end

% -------------------- Figure setup --------------------
fig = figure('Color','w');
ax = axes(fig); hold(ax,'on'); axis(ax,'equal');
grid(ax,'on'); view(ax,3);
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
title(ax,'Attitude path animation');

% Unit sphere (faint)
res = opt.sphereResolution;
[sx,sy,sz] = sphere(res);
surf(ax, sx, sy, sz, 'FaceAlpha', 0.05, 'EdgeAlpha', 0.05);

% Sun arrow
quiver3(ax, 0,0,0, sInertial(1), sInertial(2), sInertial(3), ...
    1.15, 'LineWidth', 2, 'MaxHeadSize', 0.3);
text(ax, 1.1*sInertial(1), 1.1*sInertial(2), 1.1*sInertial(3), ...
    'Sun', 'FontWeight','bold');

% Keep-out cone boundaries (draw one per distinct alpha)
alphasUnique = unique(alphaVec);
for a = alphasUnique
    ap_plotKeepOutCone(ax, sInertial, a);
end

% Plot trajectories for each sensitive axis tip
for j = 1:K
    plot3(ax, aSensIn(:,1,j), aSensIn(:,2,j), aSensIn(:,3,j), 'LineWidth', 1.2);
end

% Create animated arrows/tips for each sensitive axis
hAxis = gobjects(K,1);
hTip  = gobjects(K,1);
for j = 1:K
    hAxis(j) = quiver3(ax, 0,0,0, aSensIn(1,1,j), aSensIn(1,2,j), aSensIn(1,3,j), ...
        1.0, 'LineWidth', 2, 'MaxHeadSize', 0.25);
    hTip(j) = plot3(ax, aSensIn(1,1,j), aSensIn(1,2,j), aSensIn(1,3,j), ...
        'o', 'MarkerSize', 6, 'LineWidth', 1.5);
    text(ax, 1.05*aSensIn(1,1,j), 1.05*aSensIn(1,2,j), 1.05*aSensIn(1,3,j), ...
        sprintf('axis %d', j), 'FontSize', 9);
end

% Optional triad arrows
if opt.showTriad
    hX = quiver3(ax,0,0,0, xIn(1,1),xIn(1,2),xIn(1,3), 0.9, 'LineWidth', 1.5);
    hY = quiver3(ax,0,0,0, yIn(1,1),yIn(1,2),yIn(1,3), 0.9, 'LineWidth', 1.5);
    hZ = quiver3(ax,0,0,0, zIn(1,1),zIn(1,2),zIn(1,3), 0.9, 'LineWidth', 1.5);
end

% Limits
lim = 1.3;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax,[-lim lim]);

% -------------------- Playback indices --------------------
idxForward = 1:N;
if opt.pingpong && N >= 2
    idxForward = [1:N, (N-1):-1:2];
end

% Timing
dt = (1/opt.fps) / max(opt.speed, 1e-6);

% GIF prep
if opt.saveGif
    gifFile = opt.gifFile;
    if exist(gifFile, 'file')
        delete(gifFile);
    end
end

% -------------------- Animate --------------------
for rep = 1:max(1,opt.loops)
    for ii = idxForward
        % Update all sensitive axis arrows/tips
        angs = zeros(1,K);
        margins = zeros(1,K);

        for j = 1:K
            hAxis(j).UData = aSensIn(ii,1,j);
            hAxis(j).VData = aSensIn(ii,2,j);
            hAxis(j).WData = aSensIn(ii,3,j);

            hTip(j).XData = aSensIn(ii,1,j);
            hTip(j).YData = aSensIn(ii,2,j);
            hTip(j).ZData = aSensIn(ii,3,j);

            c = aSensIn(ii,:,j) * sInertial;   % (1x3)*(3x1)->scalar
            c = max(-1, min(1, c));
            angs(j) = acos(c);
            margins(j) = angs(j) - alphaVec(j);
        end

        % Update triad
        if opt.showTriad
            hX.UData = xIn(ii,1); hX.VData = xIn(ii,2); hX.WData = xIn(ii,3);
            hY.UData = yIn(ii,1); hY.VData = yIn(ii,2); hY.WData = yIn(ii,3);
            hZ.UData = zIn(ii,1); hZ.VData = zIn(ii,2); hZ.WData = zIn(ii,3);
        end

        % Title: show the tightest constraint
        [minMargin, jStar] = min(margins);
        title(ax, sprintf('min(angle-keepout)=%.2f deg (axis %d), angle=%.1f deg, keepout=%.1f deg', ...
            rad2deg(minMargin), jStar, rad2deg(angs(jStar)), rad2deg(alphaVec(jStar))));

        drawnow;

        if opt.saveGif
            frame = getframe(fig);
            [im, cm] = rgb2ind(frame2im(frame), 256);
            if rep == 1 && ii == idxForward(1)
                imwrite(im, cm, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', dt);
            else
                imwrite(im, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', dt);
            end
        else
            pause(dt);
        end
    end
end

end

% ========================================================================
% Helpers (prefixed to avoid naming collisions)
% ========================================================================

function A = ap_ensure3xK(A)
% Accept 3xK or Kx3; return 3xK.
if isempty(A)
    error("ASensBody cannot be empty.");
end
if size(A,1) == 3
    return;
elseif size(A,2) == 3
    A = A.';
else
    error("ASensBody must be 3xK (axes as columns) or Kx3.");
end
end

function [qDense, segId] = ap_densifyQuatPath(pathQ, samplesPerSegment)
% SLERP-interpolate each consecutive waypoint pair with fixed samples.
M = size(pathQ,1);
if M < 2
    qDense = pathQ;
    segId = 1;
    return;
end

qDense = zeros((M-1)*samplesPerSegment + 1, 4);
segId  = zeros((M-1)*samplesPerSegment + 1, 1);

row = 1;
qDense(row,:) = pathQ(1,:);
segId(row) = 1;

for i = 1:(M-1)
    q1 = pathQ(i,:);
    q2 = pathQ(i+1,:);
    ts = linspace(0,1,samplesPerSegment+1);
    ts = ts(2:end); % skip 0 (already included)
    for t = ts
        row = row + 1;
        qDense(row,:) = slerp(q1,q2,t);
        segId(row) = i;
    end
end

end

function ap_plotKeepOutCone(ax, s, alpha)
% Plot circle on the unit sphere at angular distance alpha from direction s.
s = s(:); s = s / norm(s);

% Orthonormal basis (u, v, s)
tmp = [1;0;0];
if abs(dot(tmp,s)) > 0.9
    tmp = [0;1;0];
end
u = tmp - dot(tmp,s)*s; u = u/norm(u);
v = cross(s,u);

phi = linspace(0,2*pi,250);
C = cos(alpha)*s + sin(alpha)*(u*cos(phi) + v*sin(phi));
plot3(ax, C(1,:), C(2,:), C(3,:), '--', 'LineWidth', 1.2);
end


