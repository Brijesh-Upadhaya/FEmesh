// Copyright (c) 2016 Brijesh Upadhaya

// All dimensions are in meters and rads(Pi is the only constant predefined in Gmsh)

Lc = 0.0004;            // Base char length
Z  = 0.0;               // Z-coordinate
//-------------------------------------------------------------------------------------------------//
// Stator Data 
//-------------------------------------------------------------------------------------------------//
nbs = 48;		            // No. of stator slots
dths = 2*Pi/nbs;        // Ang. shift between two slots(slot pitch)
th0s = dths/2;          // Angular pos. stator teeth
des = 0.31;             // Outer diameter of stator core
res = des/2;  	  	    // Outer radius of stator core
dis = 0.2;              // Internal diameter of stator core
ris = dis/2;		        // Internal radius of stator core
			  
hs = 0.0239;            // Slot height (H1)
h1s = 0.001;            // Height of slot opening(H11)
h3s = 0.0175;           // Height of trapezoidal portion (H13)
d1s = 0.0088;           // Diam. sup. circle (B13)
r1s = d1s/2;		        // Radius of sup. circle (R13)
h2s = hs-h1s-h3s-r1s;   // Height of section upper to slot opening (H12)
rc1s = 0.01;        	  // Radius convex rec. --> NOT used!
rc2s = 0.01;            // Radius concave rec.--> NOT used!
e1s = 0.0065;           // Slot width after slot opening (B12)
e2s = 0.0035;           // Width of slot opening (B11)

//-------------------------------------------------------------------------------------------------//
// Rotor Data (surface mounted PMSM)
//-------------------------------------------------------------------------------------------------//
nbr = 4;				        // No. of rotor magnets
dthr = 2*Pi/nbr;        // Ang. shift betw. 2 rotor magnets(Rotor slot pitch)	  
th0r = dthr/2;          // Ang. pos. of rotor gap between magnet
h1m = 0.014;            // Height of the magnet (H0)
dpm = 0.196/2;          // Rotor outer radius with magnet
der = 0.196-2*h1m;      // Rotor outer diameter: D12
rer = der/2; 		  	    // Rotor outer radius
gap = ris-rer-h1m;      // Air gap width
gap1= gap/3;            // Air-gap in three layers
dir = 0.070;          	// Internal diameter of rotor: D22
rir = dir/2; 		        // Internal radius of rotor core 
rer1 = rer+h1m+gap1;    // Inner Radius of middle layer 
rer2 = rer1 + gap1;     // Outer Radius of middle layer
beta = 0.7*dthr;        // magnet curvature

		  		  
/* 'newp' is a meta variable defining a new point number for you.
   This is mostly useful with included files. 
   'newreg'=> defines a new region number (everything that is not a point). */

pAxe = newp; Point(pAxe) = {0.0, 0.0, 0.0, 30*Lc}; // Origin
//-------------------------------------------------------------------------------------------------//
// Rotor Shaft
//-------------------------------------------------------------------------------------------------//
p1 = newp ; Point(p1) = {rir, 0.0, 0.0, 30*Lc} ;
p2 = newp ; Point(p2) = {0.0, rir, 0.0, 30*Lc} ;
Rx1 = rir * Cos(Pi/6);
Ry1 = rir * Sin(Pi/6);
Rx2 = rir * Cos(Pi/3);
Ry2 = rir * Sin(Pi/3);
pR1 = newp ; Point(pR1) = { Rx1, Ry1, 0.0, 30*Lc} ;
pR2 = newp ; Point(pR2) = { Rx2, Ry2, 0.0, 30*Lc} ;

lin1 = newreg ;	Line(lin1)   = {pAxe,p1} ;        // x-axis line in rotor
arc1 = newreg ;	Circle(arc1) = {p1,pAxe,pR1} ;    // rotorshaft outer arc1
arc2 = newreg ;	Circle(arc2) = {pR1,pAxe,pR2} ;   // rotorshaft outer arc2
arc3 = newreg ;	Circle(arc3) = {pR2,pAxe,p2} ;    // rotorshaft outer arc3
lin2 = newreg ;	Line(lin2)   = {p2, pAxe} ;       // y-axis line in rotor

// This surface includes the shaft
reg1 = newreg ; Line Loop(reg1) = {lin1,arc1,arc2,arc3,lin2};
reg2 = newreg ; Plane Surface(reg2) = {reg1}; 

//-------------------------------------------------------------------------------------------------//
// Rotor Core
//-------------------------------------------------------------------------------------------------//
p3 = newp ; Point(p3) = {rer, 0.0, 0.0, Lc} ;
p4 = newp ; Point(p4) = {0.0, rer, 0.0, Lc} ;
PMR1 = p3;
PMR2 = p4;
lin3 = newreg ; Line(lin3)   = {p1, p3} ;         // x-axis line in rotor
lin4 = newreg ;	Line(lin4)   = {p4, p2} ;         // y-axis line in rotor


// This surface includes the entire rotor core 
//reg3 = newreg ; Line Loop(reg3) = {lin3,arc2,lin4,-arc1};
//reg4 = newreg ; Plane Surface(reg4) = {reg3};
//-------------------------------------------------------------------------------------------------//
// Air-Gap (3 layers of air-gap)
//-------------------------------------------------------------------------------------------------//
// Points in the lateral side
pA1 = newp ; Point(pA1) = {rer1, 0.0, 0.0, Lc};
pA2 = newp ; Point(pA2) = {rer2, 0.0, 0.0, Lc};
pB1 = newp ; Point(pB1) = {0.0, rer1, 0.0, Lc};
pB2 = newp ; Point(pB2) = {0.0, rer2, 0.0, Lc};
p5 = newp ; Point(p5) = {ris, 0.0, 0.0, Lc} ;
p6 = newp ; Point(p6) = {0.0, ris, 0.0, Lc} ;
// Points for first arc in the middle of the Air-gap
Ax1 = rer1 * Cos(Pi/6);
Ay1 = rer1 * Sin(Pi/6);
Ax2 = rer1 * Cos(Pi/3);
Ay2 = rer1 * Sin(Pi/3);
pA11 = newp ; Point(pA11) = { Ax1, Ay1, 0.0, Lc} ;
pA12 = newp ; Point(pA12) = { Ax2, Ay2, 0.0, Lc} ;
// Points for second arc in the middle of the Air-gap
Bx1 = rer2 * Cos(Pi/6);
By1 = rer2 * Sin(Pi/6);
Bx2 = rer2 * Cos(Pi/3);
By2 = rer2 * Sin(Pi/3);
pB11 = newp ; Point(pB11) = { Bx1, By1, 0.0, Lc} ;
pB12 = newp ; Point(pB12) = { Bx2, By2, 0.0, Lc} ;
// Lateral line and Air-gap arc
lin5 = newreg ; Line(lin5) = {p3, pA1} ;
lin6 = newreg ; Line(lin6) = {pB1, p4} ;
linA6 = newreg; Line(linA6) = {pA1,pA2};
linB7 = newreg; Line(linB7) = {pB2,pB1};
linA8 = newreg; Line(linA8) = {pA2,p5} ;
linB9 = newreg; Line(linB9) = {p6,pB2} ;
// First Arc
arcA1 = newreg; Circle(arcA1) = {pA1, pAxe, pA11} ; //outer arc1
arcA2 = newreg; Circle(arcA2) = {pA11, pAxe, pA12}; //outer arc2
arcA3 = newreg; Circle(arcA3) = {pA12, pAxe, pB1} ; //outer arc3
// Second Arc
arcB1 = newreg; Circle(arcB1) = {pA2, pAxe, pB11} ; //outer arc1
arcB2 = newreg; Circle(arcB2) = {pB11, pAxe, pB12}; //outer arc2
arcB3 = newreg; Circle(arcB3) = {pB12, pAxe, pB2} ; //outer arc3

// This surface includes the Air-gap
/*reg5 = newreg ; Line Loop(reg5) = {lin5,arcA1,lin6,-arc2};*/
/*reg6 = newreg ; Plane Surface(reg6) = {reg5};*/
/*reg7 = newreg ; Line Loop(reg7) = {linA6,arcB1,linB7,-arcA1};*/
/*reg8 = newreg ; Plane Surface(reg8) = {reg7};*/
//-------------------------------------------------------------------------------------------------//
// Stator Core
//-------------------------------------------------------------------------------------------------//
p7 = newp ; Point(p7) = {res, 0.0, 0.0, 30*Lc } ;
p8 = newp ; Point(p8) = {0.0, res, 0.0 , 30*Lc } ;
Sx1 = res * Cos(Pi/6);
Sy1 = res * Sin(Pi/6);
Sx2 = res * Cos(Pi/3);
Sy2 = res * Sin(Pi/3);
pS1 = newp ; Point(pS1) = { Sx1, Sy1, 0.0, 30*Lc} ;
pS2 = newp ; Point(pS2) = { Sx2, Sy2, 0.0, 30*Lc} ;

lin7 = newreg ;	Line(lin7)   = {p5, p7} ;         // x-axis line in stator
arc1 = newreg ;	Circle(arc1) = {p7,pAxe,pS1} ;    // statoryoke outer arc1
arc2 = newreg ;	Circle(arc2) = {pS1,pAxe,pS2};    // statoryoke outer arc2
arc3 = newreg ;	Circle(arc3) = {pS2,pAxe,p8} ;    // statoryoke outer arc3
lin8 = newreg ;	Line(lin8)   = {p8, p6} ;         // y-axis line in stator

PP1 = p5 ; PPB = p6 ; 
//-------------------------------------------------------------------------------------------------//
// Rotor Slots 
//-------------------------------------------------------------------------------------------------//
D1 = dpm;
D2 = rer;
H  = h1m;     // PM Height
RPM = beta/2; // PM curvature
O2 = pAxe;    // Origin
idxR = 200;

th = th0r ;
Include "smpm_rotor.i1" ;
pMM1 = p1; 
pMM2 = p5;


// Rotor Core arc adjacent to PM
arc2 = newreg ;	Circle(arc2) = {PMR1,pAxe,pMM1} ;    // rotorcore outer arc1
arc2 = newreg ;	Circle(arc2) = {PMR2,pAxe,pMM2} ;    // rotorcore outer arc2

//-------------------------------------------------------------------------------------------------//
// Stator Slots
//-------------------------------------------------------------------------------------------------//
dth= dths; // Stator slot-pitch
D2 = ris;  // Stator inner radius
H  = hs;   // H1
R1 = r1s;  // B13/2
R2 = rc1s;
R3 = rc2s;
E1 = e1s;  // B12
E2 = e2s;  // B11
H1 = h1s;  // H11
H2 = h2s;  // H12
H3 = h3s;  // H13
O1 = pAxe; // Origin

i = 1 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP2 = p1 ; PP3 = p9 ;
i = 2 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP4 = p1 ; PP5 = p9 ;
i = 3 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP6 = p1 ; PP7 = p9 ;
i = 4 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP8 = p1 ; PP9 = p9 ;
i = 5 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP10 = p1 ; PP11 = p9 ;
i = 6 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ;
PP12 = p1 ; PP13 = p9 ;
i = 7 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP14 = p1 ; PP15 = p9 ;
i = 8 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP16 = p1 ; PP17 = p9 ;
i = 9 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP18 = p1 ; PP19 = p9 ;
i = 10 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP20 = p1 ; PP21 = p9 ;
i = 11 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP22 = p1 ; PP23 = p9 ;
i = 12 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP24 = p1 ; PP25 = p9 ;

arc1 = newreg ; Circle(arc1) = {PP1, pAxe, PP2} ;
arc1 = newreg ; Circle(arc1) = {PP3, pAxe, PP4} ;
arc1 = newreg ; Circle(arc1) = {PP5, pAxe, PP6} ;
arc1 = newreg ; Circle(arc1) = {PP7, pAxe, PP8} ;
arc1 = newreg ; Circle(arc1) = {PP9, pAxe, PP10} ;
arc1 = newreg ; Circle(arc1) = {PP11, pAxe, PP12} ;
arc1 = newreg ; Circle(arc1) = {PP13, pAxe, PP14} ;
arc1 = newreg ; Circle(arc1) = {PP15, pAxe, PP16} ;
arc1 = newreg ; Circle(arc1) = {PP17, pAxe, PP18} ;
arc1 = newreg ; Circle(arc1) = {PP19, pAxe, PP20} ;
arc1 = newreg ; Circle(arc1) = {PP21, pAxe, PP22} ;
arc1 = newreg ; Circle(arc1) = {PP23, pAxe, PP24} ;
arc1 = newreg ; Circle(arc1) = {PP25, pAxe, PPB } ;

//-------------------------------------------------------------------------------------------------//
// Domains Rotor shaft, Rotor core, Permanent Magnet, Air-gap, Stator Core
//-------------------------------------------------------------------------------------------------//
Line Loop(384) = {8, 201, -32, -31, -202, 9, -4, -3, -2};
Plane Surface(385) = {384};
Line Loop(386) = {22, 23, 24, 25, 26, -383, -364, -363, -362, -361, -360, -359, -358, -357, -382, -350, -349, -348, -347, -346, -345, -344, -343, -381, -336, -335, -334, -333, -332, -331, -330, -329, -380, -322, -321, -320, -319, -318, -317, -316, -315, -379, -308, -307, -306, -305, -304, -303, -302, -301, -378, -294, -293, -292, -291, -290, -289, -288, -287, -377, -280, -279, -278, -277, -276, -275, -274, -273, -376, -266, -265, -264, -263, -262, -261, -260, -259, -375, -252, -251, -250, -249, -248, -247, -246, -245, -374, -238, -237, -236, -235, -234, -233, -232, -231, -373, -224, -223, -222, -221, -220, -219, -218, -217, -372, -210, -209, -208, -207, -206, -205, -204, -203, -371};
Plane Surface(387) = {386};
Line Loop(388) = {201, 27, 28, 29, 30, -202, -11, -18, -17, -16, -10};
Plane Surface(389) = {388};
Line Loop(390) = {12, 19, 20, 21, 13, -18, -17, -16};
Plane Surface(391) = {390};
Line Loop(392) = {14, 371, 212, 372, 226, 373, 240, 374, 254, 375, 268, 376, 282, 377, 296, 378, 310, 379, 324, 380, 338, 381, 352, 382, 366, 383, 15, -21, -20, -19};
Plane Surface(393) = {392};
Coherence;
Show "*";
