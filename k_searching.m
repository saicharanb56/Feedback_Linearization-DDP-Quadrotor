clear; clc;

z = zeros(3);
e = eye(3);

ks = linspace(0,2,11);
K1s = [];

for a = ks
    for b = ks
        for c = ks
            for d = ks
                K = [ z e z z;
                      z z e z;
                      z z z e;
                      -a*e -b*e -c*e -d*e];
                Eig = eig(K);
                R = real(Eig);
                I = imag(Eig);
                if all(R<0)
                    K1s = [K1s; a b c d max(R) max(I)];
                end
            end
        end
    end
end

K2s = [];

for a = ks
    for b = ks
        K = [0 1; -a -b];
        Eig = eig(K);
        R = real(Eig);
        I = imag(Eig);
        if all(R<0)
            K2s = [K2s; a b max(R) max(I)];
        end
    end
end
