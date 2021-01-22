% Testing rotation procedures in CMM
% ISL 01/16/2020

% Strain displacement matrix B in lamina coords
Bl = [ -7.00388529,	           0,	           0,	  7.20040199,	           0,	           0,	-0.196516702,	           0,	           0;	
                 0,	  2.48622826,	           0,	           0,	 -2.96367567,	           0,	           0,	 0.477447403,	           0;	
        2.48622826,	 -7.00388529,	           0,	 -2.96367567,	  7.20040199,	           0,	 0.477447403,	-0.196516702,	           0;	
                 0,	           0,	 -7.00388529,	           0,	           0,	  7.20040199,	           0,	           0,	-0.196516702;	
                 0,	           0,	  2.48622826,	           0,	           0,	 -2.96367567,	           0,	           0,	 0.477447403; ];

% Elasticity tensor D
nu    = 0.5;
kappa = 0.833333;
E     = 2500000;
coef  = E / (1.0 - nu*nu);
D = coef * [ 1, nu,          0,                0,               0;
            nu,  1,          0,                0,               0;
             0,  0, 0.5*(1-nu),                0,               0;
             0,  0,          0, 0.5*kappa*(1-nu),               0;
             0,  0,          0,                0, 0.5*kappa*(1-nu) ];

% Stiffness matrix in lamina coords
Kl = Bl' * D * Bl;

% Rotation matrix
Q = [0.707889959,	0.394378312, -0.585967195;
    -0.269753948,	 0.91766973,  0.291744878;
      0.65278221, -0.0484563057,  0.755994294 ];
theta = [          Q, zeros(3, 3), zeros(3, 3);
         zeros(3, 3),           Q, zeros(3, 3);
         zeros(3, 3), zeros(3, 3),           Q ];

% Stiffness matrix in global coords
Kg = theta' * Kl * theta;

ctrl_xyz = [0.3971,  -1.4233, 9.7337;
            0.4969,  -1.2942  9.6558;
            0.4516,   1.3001, 9.8612];
ctrl_xyz = [ctrl_xyz;
            0.5 * (ctrl_xyz(1,:) + ctrl_xyz(2,:));
            0.5 * (ctrl_xyz(2,:) + ctrl_xyz(3,:));
            0.5 * (ctrl_xyz(1,:) + ctrl_xyz(3,:)) ];
