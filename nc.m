classdef nc
   properties (Constant)
        n=3;                                       %number of elements for section #3 - length 1350m
        H=1350;
        h = nc.H/nc.n;                              %length of each segment
        ad = 0.0034;                                %area of drilling collar [m^2]
        v0=10;                                     %constant speed of drilling
        %% INPUT for BIT         
        zb=5;                                       %Number of blades
        a=0.4445/2;                                 %Radius of section#3 in [m]
                                                    %Orientation of cutting face =1 at 1st iteration
        se=32e6;                                    %intrinsic specific energy (typically 32 MPa)
        e=nc.zb*nc.a*nc.se;                         %Cutting Energy       
        %% INPUT for DRILLSTRING    
        m=71000/nc.n;                               %mass of each block
        E = 200e9;                                  % Young's modulus steel
        k = nc.E*nc.ad/nc.h;                        % sprint constant (linear)
        c=535/nc.n;                                 %damping of each block
        %% Matrix
        %M-Mass matrix
        M=nc.m*eye(nc.n);                           %Matrix of drillstring mass
        %C-damping matrix
        C=nc.c*eye(nc.n);
   end
end