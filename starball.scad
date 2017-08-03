/*
 *  backy ball fractal
 *  (c) 2015. 10. 30 by SASAKI, Taroh
 */

use <maths.scad>

layer_h = .4;
extrusion_width = layer_h*1.5;
eta = 0.001;

use_lock = true;

swhex  = 1;     // for dual color printing
swpent = 1;

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
//
//  rotate 20.90516 (around y axis) degree afterwards to 
//  fit hexagon on
//  x-y plane.

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

// calculate pentagon intersection point
rat = 3/2 - sqrt(5)/2;   // ratio of petagon side to intersection distance
function penti(a,b) = a * (1-rat) + b*rat;

// 12x star (points)
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

// 12 x hex ddh around +/- z axis (0 <= x, ccw)
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

radius = 50;

//align hex1 to Z axis
hex1_center = (ddh1+ddh2+ddh3+ddh4+ddh5+ddh6)/6;
align_hex1 = quat([0,1,0], -atan2(hex1_center[0], hex1_center[2]));

// align hex2 to Z axis
// implemented as 180deg rotation of hex1 along common point
hex2_center = (och1+och2+och3+och4+och5+och6)/6;
align_hex2 = quat_mult(align_hex1, quat((hex1_center+hex2_center)/2, 180));

star_center=(ddp1+ddp2+ddp3+ddp4+ddp5)/5;
align_star = quat([0,1,0], -atan2(star_center[0], star_center[2]));

// duplicate rotating 180deg around axis
function ql_rot1d(ql, axis, angle=180) = concat(ql, [for (q = ql) quat_mult(q, quat(axis, angle))]);
// rotate twice 90 deg (to align along different axis)
function ql_rot2(ql, a1, a2, angle=90) = [for (q = ql) quat_mult(q, quat_mult(quat(a1, angle), quat(a2, angle)))];

function rot12(ql) = let(x_sym = ql_rot1d(ql_rot1d(ql, [0,1,0]),[0,0,1])) concat(x_sym, ql_rot2(x_sym, [1,0,0], [0,0,1]), ql_rot2(x_sym, [1,0,0], [0,1,0]));
// use 'mirror' rotation + 180 degree symetry
function ql_mirror(ql) = [for (q = ql) [q[0],-q[1],q[2],-q[3]]];
function rot8(ql) = ql_rot1d(ql_rot1d(concat(ql, ql_mirror(ql)), [1,0,0]), [0,0,1]);

align_p = rot12([align_star]);
align_h = concat(rot12([align_hex1]),rot8([align_hex2])); 


module align(q, rev = false) {
  multmatrix(quat_to_mat4(rev ? quat_conj(q) : q)) children();
}

// vec4_mult_mat34 is opposite to multmatrix
function align_vec(v, q, rev=false) = vec4_mult_mat34(point3_from_vec3(v), quat_to_mat4(rev ? q : quat_conj(q)));

// render element numbers + first point
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
  drawbucky();
 }

*translate([50,0,0])
  hexp(mag="left");
*hexp(mag="right");

*starp();
*translate([0,50,0])
hexp();


module stand() {
  dia=60;
  height=10;
  function pp(i) = [sin(i*360/5), cos(i*360/5)]*dia/2;
  difference() {
    eh = 20;
    translate([0,5,0]) linear_extrude(height=eh, scale=(dia/2-eh*tan(20))/(dia/2))
      polygon(points = [for(i=[0:4], j=0) j ? penti(pp(i), pp(i+2)) : pp(i)]);
    translate([0,0,radius+6]) sphere(r=radius,$fn=100);
    rotate([5,0,0]) slab(height, eh+1, dia+20);
  }
}

stand();

//prism10(ddp1, ip1, ddp2, ip2, ddp3, ip3, ddp4, ip4, ddp5, ip5);
  
module star_cut(l=[1,1,1,1,1,1]) {
  for(i = [0:2]) {
    a=[align_p[0], align_p[8],align_p[11]][i];  // stars adjacent to hex0
    sl=[[0,0,l[2],l[3],0],
	[l[0],l[1],0,0,0],
	[l[5],0,0,0,l[4]]][i];
    align(a,true)  // move to correct position
      align(align_star) drawstar(true, sl);  // star aligned to Z
  }
}

hex_h = 10;

module slab(zmin, zmax, size=200) {
  translate([-size/2,-size/2,zmin])
    cube([size,size,zmax-zmin], center=false);
}

magnet_dia=5+.5;
magnet_h = 2.1;

module line(l) {
  translate(l[0]) rotate(LineRotations(l[1]-l[0])) cylinder(r=.1,h=norm(l[1]-l[0]));
}

module hexp(mag = "none") {
  // projet to plane Z=radius
  function zproj(vl) = [for(vo=vl) let (v = align_vec(vo, align_hex1)) [v[0]/v[2], v[1]/v[2], 1]*(radius+1)];
  ohp = zproj(oh);
  ddhp = zproj(ddh);

  cutline = (mag == "left") ? [ddhp[5], ddhp[0]] : [ddhp[0], ddhp[1]];

  oh1pn = (mag == "none") ? ohp[0] : line_intersection_pp(ohp[1], ohp[0], cutline[0], cutline[1]); // extend oh1 side
  oh6pn = (mag == "none") ? ohp[5] : line_intersection_pp(ohp[4], ohp[5], cutline[0], cutline[1]); // cut oh6 side
  
  // translate(oh1pn) color("blue") cube(.2,center=true);
  // translate(oh6pn) color("red") cube(.2,center=true);
  
  ohpn = concat([oh1pn], slice(ohp, [1:4]), [oh6pn]);
  //for(i=[0:len(ohpn)-1]) translate(ohpn[i]) linear_extrude(height=.1) text(text=str(i),size=1, valign="center", halign="center");
  difference() {
    union() {
      intersection() {
	sphere(r=radius, $fn=100);
	slab(radius-hex_h, radius+1);
	union() {
	  prism_v(ddhp);
	  prism_v(ohpn);
	}
      }
    }
    locks = (mag=="none") ? [1,1,1,1,1,1] : (mag=="left") ? [0,0,1,1,1,2] : [3,1,1,1,0,0];
    align(align_hex1) star_cut(locks);
    if(mag!="none") {
      side = cutline[0]-cutline[1];
      
      v = VNORM(align_vec(ddh1, align_hex1));             // point to place magnet
      aangle = asin(v[0]*VNORM(side)[1]);                 // rotation tperpendicular point on side
      echo(v, aangle, side);
      translate(VNORM([-side[1],side[0],0]) * extrusion_width * -1)          // inset, perpendicular to side
	translate(v*(radius-magnet_dia/2-2))	                   // 2mm from surface
	  rotate(aangle, v = [side[0],side[1],0])             // align to edge
	    rotate([0,0,atan2(side[1], side[0])-90]) {
	      rotate([0,90,0]) rotate([0,0,360/16]) cylinder(d=magnet_dia/cos(180/8), h=magnet_h+.2, $fn=8);
	      translate([0, -magnet_dia/2, -10]) cube([magnet_h, magnet_dia, 10]);
	    }
    }
  }
}

star_h = 8;


//function sproj(v) = let(p=VNORM(align_vec(v, align_star))*radius,
//			h0=sqrt(radius*radius-VLENSQR([p[0], p[1], 0])))
//  p * (radius-star_h)/h0;

module starp(cut=false, locks=[1,1,1,1,1]) {
  // align and project to plane above ball ; move the plane eta down to make prism for cut slightly larger
  function zproj(vl) = [for(vo=vl) let (v = align_vec(vo, align_star)) [v[0]/v[2], v[1]/v[2], 1-(cut?eta:0)]*(radius+1)];
  // align and project to star_h plane
  function sproj(vec) = let (v = align_vec(vec, align_star)) [v[0]/v[2], v[1]/v[2], 1] * (radius-star_h);
  
  difference() {
    intersection() {
      sphere(r=radius+(cut?1:0), $fn=100);
      slab(radius-star_h, radius+1);
      prism_v(zproj([ddp1, ip1, ddp2, ip2, ddp3, ip3, ddp4, ip4, ddp5, ip5]));
    }
    if(use_lock) {
      p1 = sproj(ddp1);                                  // start tip
      v2 = VNORM(sproj(ip1) - p1);                       // unit vectors along star edges
      v3 = VNORM(sproj(ip5) - p1);
      v0 =  VNORM(v2+v3);                                // center of tip direction (pointing to star center)
      sin_a = norm(v2-v3)/2;                             // sin(tip_angle/2)
      alpha = asin(sin_a); 
      cos_a = cos(alpha);                                // cos(tip_angle/2)
      tan_a = tan(alpha);
      ofs = 2*extrusion_width;                           // star wall thickness
      lock_size = 8;                                     // size of lock triangle (centerline length)
      lock_height=layer_h*4-(cut?.3:0);                  // height of lock
      lock_cut_star = 1.5;                               // tip cut for star
      lock_cut_hex  = 2;                                 // tip cut for hexagon
      p1n = p1 + v0 * ofs / sin_a;                       // shift star tip to get ofs spacing
      p1n_hex = p1n + v0 * lock_cut_hex;                 // cut from tip on hex (cut plane intersection)
      p1n2 = p1n + v2 * lock_cut_star / cos_a;           // replace tip with line
      p1n3 = p1n + v3 * lock_cut_star / cos_a;
      p2n = p1n + v2*lock_size / cos_a;  
      p3n = p1n + v3*lock_size / cos_a;

      q1n2 = p1 + v2 * (ofs / tan_a + extrusion_width/tan_a/2);
      q1n3 = p1 + v3 * (ofs / tan_a + extrusion_width/tan_a/2);
      q2n = sproj(ip1);
      q3n = sproj(ip5);

      ofs_x = -ofs*sin_a;
      ofs_y = -ofs*cos_a;
      
      poly_base = [q1n2,q2n,q3n,q1n3];
      poly_d = [for(p=poly_base) [p[0]+ofs_x, p[1]+ofs_y*sign(p[1]), p[2]]];
      poly_top_scale = (radius-star_h+lock_height)/(radius-star_h);
      // hardcoded direction along X axis!
      poly_u = [for(p=poly_base) [p[0]*poly_top_scale+ofs_x, p[1]*poly_top_scale+ofs_y*sign(p[1]), p[2]*poly_top_scale] ];
      poly_t1 =  [for(p=poly_u) [p[0], p[1] - extrusion_width/2 * sign(p[1]), p[2]-eta]];
      poly_t2 =  [for(p=poly_u) [p[0], p[1] - extrusion_width/2 * sign(p[1]), p[2]+layer_h]];
      // norm vector of hex to the right and left
      h4n = align_vec(align_vec([0,0,1], align_h[4], true), align_star);
      h6n = align_vec(align_vec([0,0,1], align_h[6], true), align_star);
      for (i=[0:4]) {
	a = i*360/5;
	if(locks[i]) {
	  rotate([0,0,a]) {
	    difference() {
	      translate([0,0,-eta])
		intersection() {
		  union() {
		    // angled lock walls
		    polyhedron(points = concat(poly_d,poly_u),
			       faces = [[0,1,2,3], [0,4,5,1], [1,5,6,2], [2,6,7,3], [3,7,4,0], [4,7,6,5]]);
		    if(!cut)  // bridge needs something to catch on
		      polyhedron(points = concat(poly_t1,poly_t2),
				 faces = [[0,1,2,3], [0,4,5,1], [1,5,6,2], [2,6,7,3], [3,7,4,0], [4,7,6,5]]);
		    // straight walls
		    if(0) {
		      translate([0,0,radius-star_h])
			linear_extrude(lock_height) polygon(points=poly_d);
		    }
		  }
		  if(!cut)
		    translate([p1n[0]-lock_cut_star-lock_size/2,0, p1n[2]])
		      cube([lock_size, 20, 20], center=true);
		}
	      if(cut) {
		if(locks[i] == 1 || locks[i] == 2) translate(p1n_hex) multmatrix(quat_to_mat4(quat_v2([0,0,1],h4n))) slab(0,5,10);
		if(locks[i] == 1 || locks[i] == 3) translate(p1n_hex) multmatrix(quat_to_mat4(quat_v2([0,0,1],h6n))) slab(0,5,10);
	      }
	    }
	  }
	}
      }
    }
  }
}

module drawhex() {
  render()
    align(align_hex1, true) hexp(); 
}

module drawstar(cut=false, locks=[1,1,1,1,1]) {
  render()
    align(align_star, true) starp(cut, locks);
}

module
drawbucky() {
  for(a=align_p)
    align(a,true)
      color("green", .7) render() starp(); 
  for(a=align_h)
    align(a,true)
      color("blue", .7) render() hexp();
}

module prism6(p1, p2, p3, p4, p5, p6) {
  if (swhex == 1) {
    polyhedron(points = [po, p1, p2, p3, p4, p5, p6],
	       faces = [
			[1,2,3,4, 5, 6],
			[0, 2, 1], [0, 3, 2], [0, 4, 3],	
			[0, 5, 4], [0, 6, 5], [0, 1, 6]]);
  }
}

module prism_v(p) {
  sides = len(p);
  faces = concat([[for(i=[1:sides]) i]],
		 [for(i=[1:sides]) [0, (i%sides)+1, i]]);
  polyhedron(points = concat([po], p),
  	     faces = faces);
}

module prism5(p1, p2, p3, p4, p5) {
  if (swpent == 1) {
    polyhedron(points = [po, p1, p2, p3, p4, p5],
        faces = [[1, 2, 3, 4, 5],
            [0, 2, 1], [0, 3, 2], [0, 4, 3],
            [0, 5, 4], [0, 1, 5]]);
  }
}

module prism10(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
{
  if (swpent == 1) {
    polyhedron(points = [po, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 ],
	       faces = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
			[0, 2, 1], [0, 3, 2], [0, 4, 3], [0, 5, 4], [0, 6, 5], [0, 7, 6], [0, 8, 7], [0, 9, 8], [0, 10, 9],
			[0, 1, 10]]);
  }
}

