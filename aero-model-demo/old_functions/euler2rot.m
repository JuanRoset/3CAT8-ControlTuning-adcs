function R = euler2rot(roll, pitch, yaw)
    % Convert Euler angles to rotation matrix
    % Input:
    %   roll:   Rotation about x-axis (in radians)
    %   pitch:  Rotation about y-axis (in radians)
    %   yaw:    Rotation about z-axis (in radians)
    % Output:
    %   R:      3x3 rotation matrix
    
    % Precompute sin and cosine values
    cr = cos(roll);
    sr = sin(roll);
    cp = cos(pitch);
    sp = sin(pitch);
    cy = cos(yaw);
    sy = sin(yaw);
    
    % Compute rotation matrix
    R = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr;
         sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr;
           -sp,            cp*sr,            cp*cr];
end