function Rz = ROTZ(yaw)

    Rz = [cos(yaw) -sin(yaw) 0;
          sin(yaw)  cos(yaw) 0;
          0         0        1];
end