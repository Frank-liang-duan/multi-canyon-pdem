function [valid, R_ordered, D_ordered] = check_overlap(a0, R_left, R_right, D_left, D_right)
% CHECK_OVERLAP  Verify that valleys do not overlap.
%   [valid, R_ordered, D_ordered] = check_overlap(a0, R_left, R_right, D_left, D_right)
%   checks that the radius of each valley plus the radius of the previous one
%   does not exceed the prescribed spacing.
%
%   Inputs:
%       a0       - radius of central valley (fixed)
%       R_left   - vector of radii on the left side (from center outward)
%       R_right  - vector of radii on the right side
%       D_left   - spacings between successive left valleys
%       D_right  - spacings between successive right valleys
%
%   Outputs:
%       valid     - true if no overlap detected
%       R_ordered - [fliplr(R_left), R_right]
%       D_ordered - [fliplr(D_left), D_right]

valid = true;
N_side = length(R_left);
R_ordered = [];
D_ordered = [];

% Check left side
R_prev = a0;
for i = 1:N_side
    if R_prev + R_left(i) >= D_left(i)
        valid = false;
        return;
    end
    R_prev = R_left(i);
end

% Check right side
R_prev = a0;
for i = 1:N_side
    if R_prev + R_right(i) >= D_right(i)
        valid = false;
        return;
    end
    R_prev = R_right(i);
end

R_ordered = [fliplr(R_left), R_right];
D_ordered = [fliplr(D_left), D_right];
end