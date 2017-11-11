// Copyright (c) 7.11.2016 Brijesh Upadhaya

// Skeleton(cage Induction machine) provided by: Francois Henrotte

// All dimensions are in meters and rads(Pi is the only constant predefined in Gmsh)

Lc = 0.00030525;          // Base char length 
Z  = 0.0;               // Z-coordinate
//-------------------------------------------------------------------------------------------------//
// Stator Data (slot shape)
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
// Rotor Data 
//-------------------------------------------------------------------------------------------------//
nbr = 40;				        // No. of rotor slots
dthr = 2*Pi/nbr;        // Ang. shift betw. 2 rotor slots(Rotor slot pitch)	  
th0r = dths/2;          // Ang. pos. of rotor teeth
der = 0.1984;           // Rotor outer diameter: D12
rer = der/2; 		  	    // Rotor outer radius
gap = ris-rer;          // Air gap width
gap1= gap/8;            // Air-gap in three layers
dir = 0.070;          	// Internal diameter of rotor: D22
rir = dir/2; 		        // Internal radius of rotor core 
rer1 = rer + 3*gap1;    // Inner Radius of middle layer 
rer2 = rer1 + 2*gap1;   // Outer Radius of middle layer

espa = 0.0007;          // Dist. stator-pole top (H21)
dist = rer - espa;	    // Radial dist. of intersection point
d1r = 0.006;            // Diam. sup circle(B23)
r1r = d1r/2;		        // Radius of sup. cir. (R23)
d2r = 0.0025;           // Diam. inf. circle(B22=B24)
r2r = d2r/2;            // Radius of inf. cir. (R22)
d3r = 0.00585;          // Diam. middle circle(B25)
r3r = d3r/2;            // Radius of mid. cir. (R25)
hr3 = 0.0178;           // H23
hr4 = 0.0161;           // H24
hr = hr4+hr3+r2r-espa;  // Pole height (Height of rotor bar H2)
		  		  
/* 'newp' is a meta variable defining a new point number for you.
   This is mostly useful with included files. 
   'newreg'=> defines a new region number (everything that is not a point). */

pAxe = newp; Point(pAxe) = {0.0, 0.0, 0.0, 45*Lc}; // Origin
//-------------------------------------------------------------------------------------------------//
// Rotor Shaft
//-------------------------------------------------------------------------------------------------//
p1 = newp ; Point(p1) = {rir, 0.0, 0.0, 45*Lc} ;
p2 = newp ; Point(p2) = {0.0, rir, 0.0, 45*Lc} ;
Rx1 = rir * Cos(Pi/6);
Ry1 = rir * Sin(Pi/6);
Rx2 = rir * Cos(Pi/3);
Ry2 = rir * Sin(Pi/3);
pR1 = newp ; Point(pR1) = { Rx1, Ry1, 0.0, 45*Lc} ;
pR2 = newp ; Point(pR2) = { Rx2, Ry2, 0.0, 45*Lc} ;

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
// Every Pi/12 shift 
Rx33 = rer * Cos(Pi/12);
Ry33 = rer * Sin(Pi/12);
Rx3 = rer * Cos(Pi/6);
Ry3 = rer * Sin(Pi/6);
Rx44 = rer * Cos(Pi/4);
Ry44 = rer * Sin(Pi/4);
Rx4 = rer * Cos(Pi/3);
Ry4 = rer * Sin(Pi/3);
Rx55 = rer * Cos(5*Pi/12);
Ry55 = rer * Sin(5*Pi/12);
pR33 = newp ; Point(pR33) = { Rx33, Ry33, 0.0, Lc};
pR44 = newp ; Point(pR44) = { Rx44, Ry44, 0.0, Lc};
pR55 = newp ; Point(pR55) = { Rx55, Ry55, 0.0, Lc};
pR3 = newp ; Point(pR3) = { Rx3, Ry3, 0.0, Lc} ;
pR4 = newp ; Point(pR4) = { Rx4, Ry4, 0.0, Lc} ;
lin3 = newreg ; Line(lin3)   = {p1, p3} ;          // x-axis line in rotor
arc1 = newreg ;	Circle(arc1) = {p3,pAxe,pR33} ;    // rotorcore outer arc1
arc2 = newreg ;	Circle(arc2) = {pR33,pAxe,pR3} ;   // rotorcore outer arc2
arc3 = newreg ;	Circle(arc3) = {pR3,pAxe,pR44} ;   // rotorcore outer arc3
arc4 = newreg ;	Circle(arc4) = {pR44,pAxe,pR4} ;   // rotorcore outer arc4
arc5 = newreg ;	Circle(arc5) = {pR4,pAxe,pR55} ;   // rotorcore outer arc5
arc6 = newreg ;	Circle(arc6) = {pR55,pAxe,p4} ;    // rotorcore outer arc6
lin4 = newreg ;	Line(lin4)   = {p4, p2} ;          // y-axis line in rotor

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

//-------------------------------------------------------------------------------------------------//
// Stator Core
//-------------------------------------------------------------------------------------------------//
p7 = newp ; Point(p7) = {res, 0.0, 0.0, 40*Lc } ;
p8 = newp ; Point(p8) = {0.0, res, 0.0, 40*Lc } ;
Sx1 = res * Cos(Pi/6);
Sy1 = res * Sin(Pi/6);
Sx2 = res * Cos(Pi/3);
Sy2 = res * Sin(Pi/3);
pS1 = newp ; Point(pS1) = { Sx1, Sy1, 0.0, 40*Lc} ;
pS2 = newp ; Point(pS2) = { Sx2, Sy2, 0.0, 40*Lc} ;

lin7 = newreg ;	Line(lin7)   = {p5, p7} ;         // x-axis line in stator
arc1 = newreg ;	Circle(arc1) = {p7,pAxe,pS1} ;    // statoryoke outer arc1
arc2 = newreg ;	Circle(arc2) = {pS1,pAxe,pS2};    // statoryoke outer arc2
arc3 = newreg ;	Circle(arc3) = {pS2,pAxe,p8} ;    // statoryoke outer arc3
lin8 = newreg ;	Line(lin8)   = {p8, p6} ;         // y-axis line in stator

PP1 = p5 ; PPB = p6 ; 
//-------------------------------------------------------------------------------------------------//
// Rotor Slots (201,202,...)
//-------------------------------------------------------------------------------------------------//
D1 = dist; 
D2 = rer-hr4+r3r;
H  = hr;  // Slot Height
R1 = r1r; // Sup. Cir.
R2 = r2r; // Inf. Cir.
R3 = r3r; // Mid. Cir.

idxR = 200;
i = 0;
For(1:10)
  i++ ;
  th = th0r + (i - 1) * dthr;
  Include "rotor.i1" ; 
EndFor
//-------------------------------------------------------------------------------------------------//
// Stator Slots (41,42,...)
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

// Index for slots and slot-openings
idxS = 40;
idxA = 60;

i = 1 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP2 = p1 ; PP3 = p9 ; Ar1 = arc3;
i = 2 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP4 = p1 ; PP5 = p9 ; Ar2 = arc3;
i = 3 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP6 = p1 ; PP7 = p9 ; Ar3 = arc3;
i = 4 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP8 = p1 ; PP9 = p9 ; Ar4 = arc3;
i = 5 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP10 = p1 ; PP11 = p9; Ar5 = arc3;
i = 6 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ;
PP12 = p1 ; PP13 = p9; Ar6 = arc3;
i = 7 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP14 = p1 ; PP15 = p9; Ar7 = arc3;
i = 8 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP16 = p1 ; PP17 = p9; Ar8 = arc3;
i = 9 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP18 = p1 ; PP19 = p9; Ar9 = arc3;
i = 10 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP20 = p1 ; PP21 = p9; Ar10 = arc3;
i = 11 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP22 = p1 ; PP23 = p9; Ar11 = arc3;
i = 12 ; th = th0s + (i - 1) * dths ;
Include "stator.i2" ; 
PP24 = p1 ; PP25 = p9; Ar12 = arc3;

arc1 = newreg ; Circle(arc1) = {PP1, pAxe, PP2} ;
arc2 = newreg ; Circle(arc2) = {PP3, pAxe, PP4} ;
arc3 = newreg ; Circle(arc3) = {PP5, pAxe, PP6} ;
arc4 = newreg ; Circle(arc4) = {PP7, pAxe, PP8} ;
arc5 = newreg ; Circle(arc5) = {PP9, pAxe, PP10} ;
arc6 = newreg ; Circle(arc6) = {PP11, pAxe, PP12} ;
arc7 = newreg ; Circle(arc7) = {PP13, pAxe, PP14} ;
arc8 = newreg ; Circle(arc8) = {PP15, pAxe, PP16} ;
arc9 = newreg ; Circle(arc9) = {PP17, pAxe, PP18} ;
arc10 = newreg ; Circle(arc10) = {PP19, pAxe, PP20} ;
arc11 = newreg ; Circle(arc11) = {PP21, pAxe, PP22} ;
arc12 = newreg ; Circle(arc12) = {PP23, pAxe, PP24} ;
arc13 = newreg ; Circle(arc13) = {PP25, pAxe, PPB } ;

//-------------------------------------------------------------------------------------------------//
// Domains: Rotor core(8),Stator Core(12),Air-gap(9,10,11), Stator slot, Rotor slot, Rotor shaft etc
//-------------------------------------------------------------------------------------------------//
// Rotor core
Line Loop(476) = {8, 9, 10, 11, 12, 13, 14, 15, -4, -3, -2};
Plane Surface(8) = {45, 214, 227, 240, 253, 266, 279, 292, 305, 318, 476};
// Air-gap 1 (Adjacent to rotor)
Line Loop(478) = {16, 22, 23, 24, 17, -14, -13, -12, -11, -10, -9};
Plane Surface(9) = {478};
// Air-gap 2 (Middle layer)
Line Loop(479) = {18, 25, 26, 27, 19, -24, -23, -22};
Plane Surface(10) = {479};
// Air-gap 3 (Adjacent to stator)
Line Loop(480) = {20, 463, 328, 464, 340, 465, 352, 466, 364, 467, 376, 468, 388, 469, 400, 470, 412, 471, 424, 472, 436, 473, 448, 474, 460, 475, 21, -27, -26, -25};
Plane Surface(11) = {480};
// Stator core
Line Loop(481) = {28, 29, 30, 31, 32, -475, -458, -457, -456, -455, -454, -453, -452, -451, -474, -446, -445, -444, -443, -442, -441, -440, -439, -473, -434, -433, -432, -431, -430, -429, -428, -427, -472, -422, -421, -420, -419, -418, -417, -416, -415, -471, -410, -409, -408, -407, -406, -405, -404, -403, -470, -398, -397, -396, -395, -394, -393, -392, -391, -469, -386, -385, -384, -383, -382, -381, -380, -379, -468, -374, -373, -372, -371, -370, -369, -368, -367, -467, -362, -361, -360, -359, -358, -357, -356, -355, -466, -350, -349, -348, -347, -346, -345, -344, -343, -465, -338, -337, -336, -335, -334, -333, -332, -331, -464, -326, -325, -324, -323, -322, -321, -320, -319, -463};
Plane Surface(12) = {481};
// x-axis line
Physical Line(100) = {1, 8, 13, 15, 17, 25};
// Stator Yoke arc*/
Physical Line(101) = {26, 27, 28};
// Y-axis line
Physical Line(102) = {29, 18, 16, 14, 12, 5};
// first arc of the middle layer in air-gap
Physical Line(103) = {19, 20, 21};
// second arc of the middle layer in air-gap
Physical Line(104) = {22, 23, 24};
// Include the material list
Include "material.txt";
// Rotor shaft*/
Physical Surface(Steel) = {7};
// Rotor core*/
Physical Surface(Rotor_FeSi) = {8};
// Rotor Bar*/
Physical Surface(Al) = {201, 202, 203, 204, 205, 206, 207, 208, 209, 210};
// Stator Conductors*/
Physical Surface(Copper) = {41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52};
// Stator core*/
Physical Surface(Stator_FeSi) = {12};
// Air region of slot opening in stator and Layers of Air-gap
Physical Surface(Air) = { 9, 10, 11, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72};

Coherence;
Coherence Mesh;
Show "*";
