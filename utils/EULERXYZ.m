function R = EULERXYZ(angles)

    R = ROTX(angles(1))*ROTY(angles(2))*ROTZ(angles(3));

end