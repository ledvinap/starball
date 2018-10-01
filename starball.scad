/*
 *  backy ball fractal
 *  (c) 2015. 10. 30 by SASAKI, Taroh
 */

use <maths.scad>

//layer_h = .4;
//extrusion_width = layer_h*1.5;
layer_h = .2;
extrusion_width = .45;
eta = 0.001;

/*
    // ball radius
radius = 50;
hex_h = 12;
star_h = 8.2;
star_inset_radius = 1.5;             // radius delta for star inset
lock_max = 13.7;                    // lock base distance from star tip
lock_min  = 5.9;                    // lock tip distance from star tip
*/
radius = 40;
hex_h = 9.3;
star_h = 6.3;
star_inset_radius = .7;           // radius delta for star inset
lock_max = 12;                    // lock base distance from star tip
lock_min  = 4;                    // lock tip distance from star tip

// hexagonal piece
magnet_dia=5+.5;
magnet_h = 2.1;
// star piece
star_wall = 2*extrusion_width;       // thickness of lock wall on star
lock_clearance= .2;                  // clearance subtracted from hexagon lock
//star_inset_lock = extrusion_width/2; // amount of star expansion to lock it in place
star_inset_lock = .3;
    
// generate locks connecting star ray and hexagon
use_lock = true;
use_pet = true;    // allow for PET sag on bridge
mag_cover_thickness = extrusion_width * 2;

// Cartesian coordinates (from wikipedia):
//  Cartesian coordinates for the vertices of a 
//  truncated icosahedron centered at the origin are all even
//  permutations of:
//      (0, ±1, ±3ϕ)
//      (±2, ±(1+2ϕ), ±ϕ)
//      (±1, ±(2+ϕ), ±2ϕ)
//  where ϕ = (1 + √5) / 2 is the golden mean.
//  Radius squared equal to 9ϕ + 10
//  (radius ~= 4.956037)
//  The edges have length 2.

phi = (1 + sqrt(5)) / 2;

pp1 = [0, 1, 3 * phi];
pp2 = [phi, 2, 1 + 2 * phi];
pp3 = [2 * phi, 1, 2 + phi];
po = [0, 0, 0];

// original hex (before pentagon is extended to star)
oh6 = pp1;
oh5 = pp2;
oh4 = pp3;
oh3 = negy(oh4);
oh2 = negy(oh5);
oh1 = negy(oh6);

oh = [oh1,oh2,oh3,oh4,oh5,oh6];

function interp(a,b, rat) = a * (1 - rat) + b * rat;

// calculate pentagon intersection point
rat = 3/2 - sqrt(5)/2;   // ratio of petagon side to intersection distance
function penti(a,b) = interp(a,b,rat);

// star (points)
ddp1 = perm((oh1+oh6)/2);
ddp2 = perm((oh5+oh4)/2);
ddp3 = (oh4+oh5)/2;
ddp4 = negy(ddp3);
ddp5 = negy(ddp2);

// star midpoints
ip1 = penti(ddp1, ddp3);
ip2 = penti(ddp2, ddp4);
ip3 = penti(ddp3, ddp5);
ip4 = penti(ddp4, ddp1);
ip5 = penti(ddp5, ddp2);

// 12 x hex ddh1 around +/- z axis (0 <= x, ccw)
ddh6 = perm2(ip1);
ddh5 = ddp3;
ddh4 = ip3;
ddh3 = ddp4;
ddh2 = negy(ddh6);
ddh1 = perm2(ddp1);

ddh = [ddh1,ddh2,ddh3,ddh4,ddh5,ddh6];

// 8 x hex och at corner (0 <= x, 0 <= y, 0 <= z)
och1 = ddp3;
och2 = ip2;
och3 = ddp2;
och4 = perm(ip2);
och5 = perm2(ddp3);
och6 = perm2(ip2);

function negx(vec) = [-vec[0],  vec[1],  vec[2]];
function negy(vec) = [ vec[0], -vec[1],  vec[2]];
function negz(vec) = [ vec[0],  vec[1], -vec[2]];
function perm(vec) = [ vec[2], vec[0], vec[1]];
function perm2(vec) = perm(perm(vec));

//align hex1 to Z axis
hex1_center = (ddh1+ddh2+ddh3+ddh4+ddh5+ddh6)/6;
align_hex1 = quat([0,1,0], -atan2(hex1_center[0], hex1_center[2]));

// align hex2 to Z axis
// implemented as 180deg rotation of hex1 along common point
hex2_center = (och1+och2+och3+och4+och5+och6)/6;
align_hex2 = quat_mult(align_hex1, quat((hex1_center+hex2_center)/2, 180));

// align star to Z axis
star_center=(ddp1+ddp2+ddp3+ddp4+ddp5)/5;
align_star = quat([0,1,0], -atan2(star_center[0], star_center[2]));

// duplicate rotating 180deg around axis
function ql_rot1d(ql, axis, angle=180) = concat(ql, [for (q = ql) quat_mult(q, quat(axis, angle))]);
// rotate twice 90 deg (to align along different axis)
function ql_rot2(ql, a1, a2, angle=90) = [for (q = ql) quat_mult(q, quat_mult(quat(a1, angle), quat(a2, angle)))];
// generate 12 positions for stars and hex (2 per pos/neg axis)
function rot12(ql) = let(x_sym = ql_rot1d(ql_rot1d(ql, [0,1,0]),[0,0,1])) concat(x_sym, ql_rot2(x_sym, [1,0,0], [0,0,1]), ql_rot2(x_sym, [1,0,0], [0,1,0]));
// generate 8 diagonal positions for hex
// use 180 degree symetry + 'mirror' rotation
function ql_mirror(ql) = [for (q = ql) [q[0],-q[1],q[2],-q[3]]];
function rot8(ql) = ql_rot1d(ql_rot1d(concat(ql, ql_mirror(ql)), [1,0,0]), [0,0,1]);

// rotate stars to have tip0 crossing split plane (star 0 has tip0 in bottom half-sphere)
// rotation in 360/5
//       0  1  2  3  4  5  6  7  8  9 10 11
rot_p = [0, 0, 0, 0, 1, 4, 4, 1, 0, 0, 0, 0];
// rotate hex to have tip0 on equator
// rotations in 360/3
//       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
rot_h = [0, 0, 0, 0, 2, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

// stars with solid tip
//       0  1  2  3  4  5  6  7  8  9 10 11
tip_p = [1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0];

// position of hexagon magnets
//       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
mag_h = [0, 0, 0, 0, 1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 0, 0, 0, 0, 2, 1];

// hemisphere
//        0  1  2  3  4  5  6  7  8  9 10 11
hemi_p = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0];
//        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
hemi_h = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1];
// premultiply list of quaternions
function rotqlist_l(ql, rl, angle=1) = [for (i = [0:len(ql)-1]) i >= len(rl) ? ql[i] : quat_mult(quat([0,0,1], rl[i]*angle), ql[i])];

align_p = rotqlist_l(rot12([align_star]), rot_p, 360/5);
align_h = rotqlist_l(concat(rot12([align_hex1]),rot8([align_hex2])), rot_h, 360/3); 

// use quaternion to aligh children to Z axis
module align(q, reverse = false) {
  multmatrix(quat_to_mat4(reverse ? quat_conj(q) : q)) children();
}

// same as align(), but process point instead of geometry
// (vec4_mult_mat34 is opposite to multmatrix)
function align_vec(v, q, reverse=false) = vec4_mult_mat34(point3_from_vec3(v), quat_to_mat4(reverse ? q : quat_conj(q)));

function ceil_to_layer(h) = ceil(h / layer_h) * layer_h;

module square_c(size = [1, 1], center = false)
{
    if (len(center) == undef) {
	square(size, center);
    } else {
	if (len(size) != 2 || len(center) != 2) {
	    echo("Invalid square parameters", size, center);
	    #square();
	} else {
	    offset = [size[0] * center[0] / 2, size[1] * center[1] / 2];
	    translate(offset) square(size, true);
	}
    }
}

module cube_c(size = [1, 1, 1], center = false)
{
    if (len(center) == undef) {
        cube(size, center);
    } else {
        if (len(size) != 3 || len(center) != 3) {
            echo("Invalid cube parameters", size, center);
            square();
        } else {
            offset = [size[0] * center[0] / 2, size[1] * center[1] / 2, size[2] * center[2] / 2];
            translate(offset) cube(size, true);
        }
    }
}


// render full ball with element numbers (for each shape) + mark first vertex (others are clockwise)
*union() { 
  for(i = [0:len(align_h)-1])
    align(align_h[i],true) {
      translate([0,0,radius*1.05]) linear_extrude(.1) text(text=str(i,"."), font = "Liberation Sans", size=radius/10, valign="center", halign="center");
      mark = align_vec(ddh1, align_hex1);
      translate(VNORM([mark[0]*.8, mark[1]*.8, mark[2]])*radius*1.05) cube(1, center=true);
    }

  for(i = [0:len(align_p)-1])
    align(align_p[i],true) {
      translate([0,0,radius+3]) linear_extrude(.1) text(text=str(i,"."), font = "Liberation Sans", size=radius/10, valign="center", halign="center");
      mark = align_vec(ddp1, align_star);
      translate(VNORM([mark[0]*.8, mark[1]*.8, mark[2]])*radius*1.05) cube(1, center=true);
    }
  intersection() {
      drawball(half=1);
      //      translate([-1,0,0])cube(120);
  }
  
}

*union() {
  translate([0,50,0]) drawhex(mag="right");
  drawhex();
  translate([0,-50,0]) drawhex(mag="left");
 }

*union() {

    drawstar();
    
    translate([0,55,0]) drawstar(locks=[0,1,1,1,1]); 
 }


stand();



module stand() {
  dia=56;
  height=10;
  function pp(i) = [sin(i*360/5), cos(i*360/5)]*dia/2;
  difference() {
    eh = 20;
    translate([0,5,0]) linear_extrude(height=eh, scale=(dia/2-eh*tan(20))/(dia/2))
      polygon(points = [for(i=[0:4], j=0) j ? penti(pp(i), pp(i+2)) : pp(i)]);
    translate([0,0,radius+6]) sphere(r=radius,$fn=100);
    rotate([5,0,0]) slab_z(height, eh+1, dia+20);
  }
}

// align stars adjacent to hex 0, enable adjacent locks
module star_cut(locks = [1, 1, 1, 1, 1, 1]) {
  l = locks;
  for(i = [0:2]) {
    a=[align_p[0], align_p[8],align_p[11]][i];  // stars adjacent to hex0
    sl=[[0,    0,    l[2], l[3], 0    ],
	[l[0], l[1], 0,    0,    0    ],
	[l[5], 0,    0,    0,    l[4]]][i];
    align(a, reverse=true)  // move to correct position
      starp(cut=true, locks=sl);  // star aligned to Z
  }
}

// hexa piece, magnet can be none/left/right
module hexp(mag = "none") {
    // align to Z and project to plane Z=radius+1
    function zproj(vl) = [for(vo=vl) let (v = align_vec(vo, align_hex1)) v / v.z * (radius + 1)];
    ohp = zproj(oh);
    ddhp = zproj(ddh);

    // line between hexagons, used to cut/extend edges (unused without mag)
    cutline = (mag == "left") ? [ddhp[5], ddhp[0]] : [ddhp[0], ddhp[1]];

    // extend odd hexagon sides (0,2,4)
    function hex_extend(vl) = [for(i=[0:len(vl)-1]) let(l=len(vl), ip=(i % 2) ? modp(i + 1, l) : modp(i - 1, l)) interp(vl[ip], vl[i], 1.25)];
    ohp_1 = hex_extend(ohp);
  
    oh1pn = (mag == "none") ? ohp_1[0] : line_intersection_pp(ohp_1[1], ohp_1[0], cutline[0], cutline[1]); // oh1 side
    oh6pn = (mag == "none") ? ohp_1[5] : line_intersection_pp(ohp_1[4], ohp_1[5], cutline[0], cutline[1]); // oh6 side
  
    ohpn = concat([oh1pn], slice(ohp_1, [1:4]), [oh6pn]);
    // for(i=[0:len(ohpn)-1]) translate(ohpn[i]) linear_extrude(height=.1) text(text=str(i),size=1, valign="center", halign="center");

    difference() {
	// basic shape
	union() {
	    intersection() {
		sphere(r=radius, $fn=100);
		slab_z(radius-hex_h, radius+1);
		pyramid(ohpn);
	    }
	}
	// cut edges, with correct locks
	locks = (mag=="none") ? [1,1,1,1,1,1] : (mag=="left") ? [0,0,1,1,1,2] : [3,1,1,1,0,0];
	align(align_hex1)
	    star_cut(locks);
	// hole for magnet
	if(mag != "none") {
	    side = cutline[0]-cutline[1];
      
	    v = VNORM(ddhp[0]);                            // point to place magnet
	    translate(VNORM([-side[1],side[0],0]) * mag_cover_thickness * -1)      // inset, perpendicular to side in XY, single perimeter
		translate(v*(radius-magnet_dia/2-2)) {	                           // 2mm from sphere surface
		    aangle = asin(v.x*VNORM(side).y);                            // angle of edge in place of magnet
		    rotate(aangle, v = [side.x,side.y,0])                        // align to edge
			rotate([0,0,atan2(side.y, side.x)-90]) {                   // align to side in XY
			    // octagon for magnet, with some reserve in height
			    rotate([0,90,0]) rotate([0,0,360/16]) cylinder(d=magnet_dia/cos(180/8), h=magnet_h+.2, $fn=8);
			    // access hole
			    cube_c([magnet_h, magnet_dia, 10], center = [1, 0, -1]);
			}
		}
	}
    }
}

// star piece
module starp(cut=false, locks=[1,1,1,1,1]) {
    // align and project to plane above ball ; move the plane eta down in Z to make pyramid for cut slightly larger
    function zproj(vl) = [for(vo=vl) let (v = align_vec(vo, align_star)) (v * ((radius + 1) / v.z)) - [0,0,cut ? eta:0]];
    // align and project to star_h plane / cut clearance lower
    star_z_ref = radius - star_h;
    star_z = star_z_ref - (cut ? lock_clearance : 0);
    function sproj(vec) = let (v = align_vec(vec, align_star)) v * (star_z / v.z);

    // offset point list by v toward x
    function ofs_poly_x(pl, both, left, right) = let(l = left != undef ? left : negy(both),
                                                     r = right != undef ? right : both)
	[for (p=pl) p.y < 0 ? [p.x+l.x, p.y+l.y, p.z+l.z] : [p.x+r.x, p.y+r.y, p.z+r.z]];

    // move point by ofs perpendicular to origin (pentagon tip)
    function pentexpand_p(p, ofs) = let(v1 = VNORM([p.x, p.y]), o = [v1.y, -v1.x, 0] * ofs) p + o;
    // move pentagon point away from origin
    function pentexpand_m(p, ofs) = let(o = VNORM([p.x, p.y, 0]) * ofs) p + o;
    function pentexpand(p, ofs) =
	[for(i=[0:4], d=[-1,0,1])
	    let(o = (locks[i] && locks[modp(i+1,5)]) ? ofs : eta)   // no inset around lock
		d ? pentexpand_p(p[d<0 ? i*2 : modp(i*2+2, 10)], o*d/cos(18)) :
		    pentexpand_m(p[i*2+1], o/sin(108/2))]; // left, mid, right
    difference() {
	intersection() {
	    if(!cut) sphere(r=radius+(cut?1:0), $fn=100);
	    slab_z(star_z, radius+2);
	    pent = zproj([ddp1, ip1, ddp2, ip2, ddp3, ip3, ddp4, ip4, ddp5, ip5]);
	    union() {
		pyramid(pent);
		if(star_inset_lock > 0) {
		    p_u = pentexpand(pent, star_inset_lock);
		    p_d = pentexpand(pent*.5, star_inset_lock);  // somewhere under star base
		    intersection() {
			prism(base = p_d, top = p_u);
			sphere(r=radius-star_inset_radius, $fn=100);
		    }
		}
	    }
	}
	if(use_lock) {
	    for (i=[0:4]) if (locks[i]) {
		// lock symmetrical around X axis
		p1 = sproj(ddp1);                                  // star tip
		v2 = VNORM(sproj(ip1) - p1);                       // unit vector along star edge
		v0 = [-1,0,0];                                     // center of tip direction (pointing to star center)
		sin_a = v2.y;                                      // sin(tip_angle/2)
		alpha = asin(sin_a);                               // tip halfangle
		cos_a = cos(alpha);                                // cos(tip_angle/2)
		tan_a = tan(alpha);
	    
		inset_r = locks[modp(i+1,5)] && locks[i];
		inset_l = locks[modp(i-1,5)] && locks[i];
		ofs_r = star_wall - (inset_r ? star_inset_lock : 0);                 // star wall thickness (inset is added from outside)
		ofs_l = star_wall - (inset_l ? star_inset_lock : 0);                 // star wall thickness (inset is added from outside)
		maxofs = max(ofs_l, ofs_r);
		lock_max_a = lock_max - (cut ? .4 : 0);     // lock base distance from star tip
		lock_min_a = lock_min + (cut ? .4 : 0);      // lock tip distance from star tip
		lock_height = cut ? ceil_to_layer(1.4) : ceil_to_layer(2); // height of lock (smaller when cutting hex)
		if(lock_min_a < min(ofs_l,ofs_r) / sin_a) {
		    echo("lock_min_a is too small", lock_min = lock_min_a, min =  min(ofs_l,ofs_r) / sin_a);
		}
		    

		q1n2 = p1 + v2 * (maxofs + extrusion_width + eta) / tan_a;  // replace tip with line (tip must not degenerate after aplying bridge ofs
		q2n = sproj(ip1);                                  // point far enough on tip edge

		ofsv_r = [0, -cos_a, 0] * ofs_r;
		ofsv_l = [0,  cos_a, 0] * ofs_l;
		
		poly_base = [q1n2, q2n, negy(q2n), negy(q1n2)];                                  // polygon matching star base, with tip cut
		// offset point list by ofs toward X axis
		poly_d = ofs_poly_x(poly_base*(1-eta/star_z), left = ofsv_l, right = ofsv_r);    // bottom of prism, shift down slightly
		poly_top_scale = (star_z_ref + lock_height) / star_z;                            // scale for top of prism
		poly_u = ofs_poly_x(poly_base*poly_top_scale, left = ofsv_l, right = ofsv_r);    // scale top of prism then apply offset
		poly_t1 = ofs_poly_x(poly_u, both = [0, -extrusion_width/2, 0]);                 // step to make bridging over lock easier
		poly_t2 = ofs_poly_x(poly_u, both = [0, -extrusion_width, layer_h]);             // second step for PET
	    
		// norm vector of hex to the right and left
		h4n = align_vec(align_vec([0,0,1], align_h[4], true), align_star);
		h6n = align_vec(align_vec([0,0,1], align_h[6], true), align_star);

		rotate([0,0,i*360/5]) {                  // star tip direction
		    difference() {
			intersection() {
			    union() {
				// angled lock walls
				//echo(i=i, poly_d=poly_t1,poly_u=poly_t2);
				prism(base = poly_d, top = poly_u);
				if(!cut) {  // bridge needs something to catch on
				    prism(base = poly_t1, height=layer_h);
				    if(use_pet && layer_h < .4)
					prism(base = poly_t2, height=layer_h);
				}
			    }
			    union() {
				// cut star lock to correct size (perpendicular to star bed)
				// trim hex lock if necessary (single-sided lock may interfere)
				d = lock_max_a-lock_min_a;
				translate(p1 + [-lock_max_a, 0, 0])
				    slab_x(0, d/2, 10);
				translate(p1 + [-lock_min_a, 0, 0]) rotate([0,-30,0])
				    slab_x(-d/2, 0, 10);
			    }
			}
			if(cut) {
			    // cut parallel to hex bed
			    if(locks[i] == 1 || locks[i] == 2)
				translate(p1 + v0 * lock_min_a) multmatrix(quat_to_mat4(quat_v2([0,0,1],h4n))) { slab_z(0,5,20); }
			    if(locks[i] == 1 || locks[i] == 3)
				translate(p1 + v0 * lock_min_a) multmatrix(quat_to_mat4(quat_v2([0,0,1],h6n))) { slab_z(0,5,20); }
			    // cut perpendicular hex bed (use XZ projection of one adjacent hex)
			    translate(p1 + v0 * lock_max_a) rotate([0,atan2(h4n.x,h4n.z),0]) slab_x(-10,0,20);
			}
		    }
		}
	    }
	}
    }
}

module drawhex(mag="none") {
  render(15)
    hexp(mag); 
}

module drawstar(locks=[1,1,1,1,1]) {
  render(15)
    starp(cut=false, locks=locks);
}

module drawball(star, hex, half=0) {
    starl = star == undef ? [0:len(align_p)-1] : star;
    hexl = hex == undef ? [0:len(align_h)-1] : hex;
    for(i = starl) {
	if(half == 0 || !!hemi_p[i] == (half > 0)) {
	    locks = tip_p[i] ? [0,1,1,1,1] : [1,1,1,1,1];
	    align(align_p[i],true)
		color("green", .7) render(15) starp(locks = locks);
	}
    }
    for(i = hexl) {
	if(half == 0 || !!hemi_h[i] == (half > 0)) {
	    mag = ["none","left","right"][mag_h[i]];
	    align(align_h[i],true)
		color("blue", .9) render(10) hexp(mag);
	}
    }
}

// pyramid with tip on origin, list od points forms base
// base should be planar
module pyramid(p) {
  sides = len(p);
  faces = concat([[for(i=[1:sides]) i]],
		 [for(i=[1:sides]) [0, i % sides + 1, i]]);
  polyhedron(points = concat([po], p),
  	     faces = faces);
}

// prism starting at p[0].z
module prism(base, height, top) {
    if(top != undef) {
	if(len(base) != len(top)) {
	    echo("Base must match top", base = base, top = top);
	    #cube();
	} else {
	    l = len(base);
	    polyhedron(points = concat(base, top),
		       faces = concat([listr([0:l-1])], // base
				      [listr([2*l-1:-1:l])], // top
				      [for(i=[0:l-1]) [modp(i+1,l), modp(i,l), modp(i,l)+l, modp(i+1,l)+l]])); // sides
	}
    } else {
	ofs = len(height) == undef ? [0,0,height] : height;
	prism(base = base, top = vvec_plus(base, ofs));
    }
}


module slab_z(zmin, zmax, size=200) {
  translate([-size/2,-size/2,zmin])
    cube([size, size, zmax - zmin], center=false);
}

module slab_x(xmin, xmax, size=200) {
  translate([xmin, -size/2,-size/2])
    cube([xmax-xmin, size, size], center=false);
}



module line(l) {
  translate(l[0]) rotate(LineRotations(l[1]-l[0])) cylinder(r=.1,h=norm(l[1]-l[0]));
}