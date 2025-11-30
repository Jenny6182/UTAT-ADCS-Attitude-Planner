function U = potential(q, qf, refs, obs)
% Calculate total potential
arguments (Input)
    q
    qf
    refs
    obs
end

arguments (Output)
    U
end

exc_cam_h = deg2rad(15);
exc_str_h = deg2rad(40);
sun_width_h = deg2rad(0.25);
moon_width_h = deg2rad(0.25);
ring = deg2rad(5);
a = [qrot(q, refs(1,:))',
     qrot(q, refs(2,:))',
     qrot(q, refs(3,:))']';

f = [attractor(5.0, q, qf);
     repulsor_quadratic(a(1,:), obs(1,:), exc_cam_h + sun_width_h + ring, 100.0);
     repulsor_quadratic(a(2,:), obs(2,:), exc_str_h + sun_width_h + ring, 100.0);
     %repulsor_quadratic(a(3,:), obs(3,:), exc_str_h + moon_width_h, 1.0)];
     repulsor_gaussian(a(3,:), obs(3,:), 1.0, 0.2)];
U = sum(f);
end