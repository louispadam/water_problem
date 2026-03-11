function return_data = angle_to_perp(theta)
%ANGLE_TO_PERP Given an angle (in radians) return a unit vector
%corresponding to a perpendicular direction (hence theta + pi/2).
%
%last updated 03/09/26 by Adam Petrucci

    return_data = angle_to_vec(theta + pi/2);
    
end