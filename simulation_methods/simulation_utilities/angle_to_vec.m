function return_data = angle_to_vec(theta)
%ANGLE_TO_VEC Given an angle (in radians) return a unit vector pointing in
%that direction.
%
%last updated 03/09/26 by Adam Petrucci

    return_data = [cos(theta); sin(theta)];

end