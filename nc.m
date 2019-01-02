classdef nc
   properties (Constant) 
        n=3;                                        % number of elements for section #3 - length 1350m
        E = 200e9;                                  % Young's modulus steel
        H=1350;
        h = nc.H/nc.n;                              % length of each segment
        v0=100;                                      % constant speed of drilling
        B=0.825;                                    % bouyancy factor
        %% INPUT for BIT         
        zb=5;                                       % Number of blades
        a=0.4445/2;                                 % Radius of section#3 in [m]
                                                    % Orientation of cutting face =1 at 1st iteration
        se=32e6;                                    % intrinsic specific energy (typically 32 MPa)
        e=nc.zb*nc.a*nc.se;                         % Cutting Energy       
        %% INPUT for DRILLSTRING    
        m=71000/nc.n;                               % mass of each block        M*/n
        k = nc.n*524287                             % spring constant           n*K
        c=535/nc.n;                                 % damping of each block     C/n
        %% Matrix
        %M-Mass matrix
        M=nc.m*eye(nc.n);                           % Matrix of drillstring mass
        %C-damping matrix
        C=nc.c*eye(nc.n);
   end
end