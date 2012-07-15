.arm
.arch armv7-a
.fpu neon
.syntax unified
.global main

# previously known as intro_do, main loop
main:
	blx main_thumb

main_thumb:
.thumb
	# init GFX
	mov r1, $0
	movt r1, $0x1002

	movw r3, $0x3F9C
	movt r3, $0x3F1F
	str r3, [r1, $0x0]

	movw r3, $0x61DF
	movt r3, $0x090B
	str r3, [r1, $0x4]

	mov r3, $0x1800
	movt r3, $0x067F
	str r3, [r1, $0x8]

	mov r2, $0
	movt r2, $0x6002 
	str r2, [r1, $0x10]

	movw r3, $0x082B
	str r3, [r1, $0x18]

	# init FPU
	ldr r1, =0x40000000	@ VFPEnable
	fmxr fpexc, r1

	mov r0, $0
	mov r3, $0x12c000


	bl	main_loop

	.redbars:
		mov r1, r0
		and r1, r1, $0xFF 
		str r1, [r2, r0]
		add r0, r0, $4
		cmp r0, r3
	bne .redbars

	#yloop:
	#	xloop:
	#		steploop:
	#		bne steploop
	#	bne xloop
	#bne yloop

	# HALT, hammerzeit!
	lock:	
		b lock

# TODO implement me
init_view_port:
bx lr

# args x (s0), a (s1) , b (s2)
# return (x<a)?a:(x>b)?b:x;
clamp:
	bl max
	vmov s2, s1
	bl min
	bx lr

# args: x (s0), a (s1)
# return (x<a)?x:a;
min:
	# if (x < a) return x;
	vcmp.f32 s0, s1
	blt end_min

	# else return a;
	vmov s0, s1
end_min:
	bx lr

# args: x (s0), a (s1)
# return (x>a)?x:a;
max:
	# if (x > a) return x;
	vcmp.f32 s0, s1
	bgt end_max

	# else return a;
	vmov s0, s1
end_max:
	bx lr

# args: x (s0), y (s1)
# return sqrtf(x*x+y*y);
length2:
	vmul.f32 s0, s0
	vmul.f32 s1, s1
	vadd.f32 s0, s1
	vsqrt.f32 s0, s0
bx lr

# args: x (s0), y (s1), z (s2)
# return sqrtf(x*x+y*y+z*z);
length3:
	vmul.f32 s0, s0
	vmul.f32 s1, s1
	vmul.f32 s2, s2
	vadd.f32 s0, s1
	vadd.f32 s0, s2
	vsqrt.f32 s0, s0
bx lr

# args: s0,s1,s2
normalize:
	vpush.f32 {s0,s1,s2}
	bl length3

	vcvt.f64.f32 d16,s0
	vrecpe.f32 d16, d16
	vcvt.f32.f64 s3,d16
	
	vpop.f32 {s0,s1,s2}

	vmul.f32 s0, s3
	vmul.f32 s1, s3
	vmul.f32 s2, s3

	bx lr

# args: posx, posy, posz, r
# return length(posx,posy,posz) - r;
spheredist:
	bl length3
	vsub.f32 s0, s3
bx lr

mix:
bx lr

# args: s0=px, s1=py, s2=pz, s3=bx, s4=by, s5=bz, s6=r
udroundbox:
	# preserve some args
	vmov s8, s1
	vmov s9, s2

	# arg1 to max is always 0
	vsub.f32 s1, s1
	
	# px = max(abs(px)-bx, 0.0)
	vabs.f32 s0, s0	
	vsub.f32 s0, s3
	bl max
	vmov s7, s0

	# py = max(abs(py)-by, 0.0)
	vabs.f32 s0, s8
	vsub.f32 s0, s4
	bl max
	vmov s8, s0

	# pz = max(abs(pz)-bz, 0.0)
	vabs.f32 s0, s9
	vsub.f32 s0, s5
	bl max
	
	# t = length3(px, py, pz)
	vmov s2, s0 
	vmov s0, s7
	vmov s1, s8
	bl length3

	# return t-r
	vsub.f32 s0, s6
bx lr

# args: px, px, pz, tx, ty
sdtorus:
	bl length2
	vmov s5, s0
	#tmp=bla
	vsub.f32 s5, s3
	# tmp-=tx
	vmov s0, s5
	# arg0=bla-tx
	vmov s1, s2
	# arg1=pz
	bl length2
	vsub.f32 s0, s4
bx lr

# args: px, py, pz
dist:
	# preserve args
	vmov s10, s0
	vmov s11, s1
	vmov s12, s2

	# res1 = sdtorus(px-torus[0], py-torus[1], pz-torus[2], 3.0, 1.0)
	ldr r0, =torus
	vldr.f32 s3, [r0]
	vldr.f32 s4, [r0,#4]
	vldr.f32 s5, [r0,#8]

	# px -= torus[0]
	vsub.f32 s0, s3
	# py -= torus[1]
	vsub.f32 s1, s4
	# pz -= torus[2]
	vsub.f32 s2, s5

	vldr.f32 s3, [r0,#12]
	vldr.f32 s4, [r0,#16]
	bl sdtorus

	# preserve result
	vmov s13, s0

	# res2 = udroundbox(px-box[0], py-box[1], pz-box[2], 0.75, 3.0, 0.5, 1.0)
	vmov s0, s10
	vmov s1, s11
	vmov s2, s12

	ldr r1, =box
	vldr.f32 s3, [r0]
	vldr.f32 s4, [r0,#4]
	vldr.f32 s5, [r0,#8]

	# px -= box[0]
	vsub.f32 s0, s3
	# py -= box[1]
	vsub.f32 s1, s4
	# pz -= box[2]
	vsub.f32 s2, s5

	vldr.f32 s3, [r0,#12]
	vldr.f32 s4, [r0,#16]
	vldr.f32 s5, [r0,#20]
	vldr.f32 s6, [r0,#24]

	bl udroundbox

	# return min(res1, res2)
	vmov s1, s0
	vmov s0, s13
	bl min
bx lr

# args: x1,x2,y1,y2,z1,z2
# return x1*x2 + y1*y2 + z1*z2
dot:
	vmul.f32 s0, s1
	vmul.f32 s2, s3
	vmul.f32 s4, s5

	vadd.f32 s0, s2
	vadd.f32 s0, s4
bx lr

# args: s0, s1, s2, s3
normal_at_point:
	VPUSH.F32	{s4,s5,s6,s7,s8,s9,s10,s11,s12}

	VMOV.F32	s4,s0
	VMOV.F32	s5,s0
	VMOV.F32	s6,s0

	VMOV.F32	s7,s1
	VMOV.F32	s8,s1
	VMOV.F32	s9,s1

	VMOV.F32	s10,s2
	VMOV.F32	s11,s2
	VMOV.F32	s12,s2

	@ load epsilon
	LDR r0,=epsilon
	VLDR.F32 s0,[r0]
	
	@ += epsilon
	VADD.F32	s4,s0
	VADD.F32	s8,s0
	VADD.F32	s12,s0

	@ preserve s3
	VPUSH.F32	{s3}

	@ central diff calc
	VMOV.F32	s0,s4
	VMOV.F32	s1,s5
	VMOV.F32	s2,s6

	BL	dist
	@ preserve result
	VPUSH.F32	{s0}

	VMOV.F32	s0,s7
	VMOV.F32	s1,s8
	VMOV.F32	s2,s9

	BL	dist
	@ preserve result
	VPUSH.F32	{s0}
	
	VMOV.F32	s0,s10
	VMOV.F32	s1,s11
	VMOV.F32	s2,s12
	
	BL	dist
	VMOV.F32	s2,s0
	VPOP.F32	{s1}
	VPOP.F32	{s0}

	@ subtract d
	VPOP.F32	{s3}
	VSUB.F32	s0,s3
	VSUB.F32	s1,s3
	VSUB.F32	s2,s3

	@ normalize

	BL	normalize
bx lr


main_loop:

	LDR r0,=floats
	VLDR.32 s31,[r0] 	@ dir.x
	VLDR.32 s30,[r0] 	@ dir.y
	VLDR.32 s29,[r0]	@ dir.z

	MOV r12, $0		@ i
	MOV r11, $0		@ j

	VLDR.32 s28,[r0]	@ step

				@ s27 = ray.x, s26 = ray.y, s25 = ray.z

	VLDR.32 s24,[r0]	@ color.r
	VLDR.32 s23,[r0]	@ color.g
	VLDR.32 s22,[r0]	@ color.b

				@ s21 = specular

	MOV r10, $0		@ colorArray.r
	MOV r9, $0		@ colorArray.g
	MOV r8, $0		@ colorArray.b

	@ s20 = d
	LDR r1, =viewport
	VLDR.32 s19,[r1]	@ dx
	VLDR.32 s18,[r1,#4]	@ dy

	@ s17, s16, s15 = N
	@ s14, s13, s12 = L
	@ s11, s10, s9  = H
	

outerloop:

	@ increment i
	ADD	r12,#1
	
	@ j = 0
	MOV	r11,$0

innerloop:

	@ dir[0] = (j*dx) - 1.0;
	VMOV.F32	r11,s0			@ transfer j to fp register
	VCVT.F32.U32	s0,s0			@ convert j to floating point
	VMOV.F32	s31,#-1.0		@ s31 = -1.0
	VMLA.F32	s31,s0,s19		@ s31 = j*dx - 1.0
	
        @ dir[1] = (i*dy) - 1.0
	VMOV.F32	r12,s0			@ transfer i to fp register
	VCVT.F32.U32	s0,s0			@ convert i to floating point
	VMOV.F32	s30, #-1.0		@ s30 = -1.0
	VMLA.F32	s30,s0,s18		@ s30 = i*dy - 1.0

	@ dir[2] = -eye[2];
	VMOV.F32	s29, #-1.0		@ s29 = -1.0

	@ ray.xyz = eye.xyz
	VSUB.F32	s27, s27		@ ray[0] = 0.0
	VSUB.F32	s26, s26		@ ray[1] = 0.0
	VMOV.F32	s25, #-1.0		@ ray[2] = -1.0

	@ step = 0.0
	VSUB.F32	s28,s28			@ step = 0.0
	
raymarch:
	@ d = (dist(ray[0], ray[1], ray[2]));
	VMOV.F32	s0,s27			@ s0 = ray[0]
	VMOV.F32	s1,s26			@ s1 = ray[1]
	VMOV.F32	s2,s25			@ s2 = ray[2]

	BL dist					@ call dist()
	VMOV.F32	s20,s0			@ d = result

	@ ray hit
	VCMP.F32	s20,#0			@ d > 0 ?
	BGT		nohit			@ skip hit code

hit:
	@ NormalAtPoint(ray[0], ray[1], ray[2], d, &N[0]);
        VMOV.F32        s0,s27                  @ s0 = ray[0]
        VMOV.F32        s1,s26                  @ s1 = ray[1]
        VMOV.F32        s2,s25                  @ s2 = ray[2]
	VMOV.F32	s3,s20			@ s3 = d
	BL		normal_at_point
	
	@ L[0] = -ray[0]; L[1] = -ray[1]; L[2] = -ray[2] - 10.0
	VMOV.F32	s14,#-1.0		@ L.x = -1.0
	VMOV.F32	s12,#-10.0		@ L.z = -10.0

	VMUL.F32	s1,s26,s14		@ L.y = ray[1] * -1.0
	VMLA.F32	s2,s25,s14		@ L.z = -10.0 + (ray[2] * -1.0)
	VMUL.F32	s0,s27,s14		@ L.x = -1.0 * ray[0]

	@ Normalize L
	BL normalize
	
	@ H = L - ray
	@ VSUB.F32	s0,s27
	@ VSUB.F32	s1,s26
	@ VSUB.F32	s2,s25	

	@ Normalize H
	@ BL normalize

	@ calculate diffuse term
	VMOV.F32	s3,s17		@ s3 = N.x
	VMOV.F32	s4,s16		@ s4 = N.y
	VMOV.F32	s5,s15		@ s5 = N.z

	BL 	dot

	VSUB.F32	s1,s1		@ s1 = 0.0
	VMOV.F32	s2,#1.0		@ s2 = 1.0

	BL 	clamp

	@ color = colorDiffuse * diffuseTerm

	LDR r0,=diffuse
	
        VLDR.F32 s3,[r0]        @ diffuse.r
        VLDR.F32 s4,[r0,#4]     @ diffuse.g
        VLDR.F32 s5,[r0,#8]     @ diffuse.b

	VMUL.F32 s3,s0	
	VMUL.F32 s4,s0
	VMUL.F32 s5,s0

	@ specular = dot(N,H)
	@ specular = pow(specular, 50)

	@ color += specular*(1.0 - color)
	@ colorArray = color * 255

	VMOV.F32 s0,#16.0
	VMUL.F32 s0,s0

	VMUL.F32 s3,s0	
	VMUL.F32 s4,s0
	VMUL.F32 s5,s0
		
	VCVT.U32.F32 s3,s3
	VCVT.U32.F32 s4,s4
	VCVT.U32.F32 s5,s5

	VMOV	r10,s3
	VMOV	r9,s4
	VMOV	r8,s5

	BL doneraymarch

nohit:
	@ step += d
	VADD.F32	s28,s20			@ step += d

	@ colorArray.rgb = 0
        MOV r10, $0     		        @ colorArray.r
        MOV r9, $0	             		@ colorArray.g
        MOV r8, $0              		@ colorArray.b

	@ ray = step * dir
	VMUL.F32	s27,s31,s28		@ ray[0] = step * dir[0]
	VMUL.F32	s26,s30,s28		@ ray[1] = step * dir[1]
	VMUL.F32	s25,s29,s28		@ ray[2] = step * dir[2]

	@ loop again
	BL 		raymarch

doneraymarch:	
	
	@ plot pixel into buffer
        MOVT 	r2, $0x6002 
	MOV	r3, #640
	MUL	r3,r12,r3
	ADD	r3,r11
	
	@ alpha
	MOV	r4, $0xFF
	STR	r4, [r2, r3]

	@ R
	ADD	r3, $0x4
	STR	r10, [r2, r3]
	@ G
	ADD	r3, $0x4
	STR	r3, [r2, r3]
	@ B
	ADD	r3, $0x4
	STR	r8, [r2, r3]
	

doneinnerloop:
	MOV	r7,#640
	CMP	r11,r7
	BGT	outerloop

	@ increment j
	ADD	r11,#1
	BL	innerloop

doneouterloop:
	MOV	r7,#480
	CMP	r12,r7
	BLT	outerloop


@ LOOP DONE
	B lock



# last 2 are torus() args
torus:
	.float 4.0, 0.0, 9.0, 3.0, 1.0

box:
	.float -3.0, 0.0, 10.0, 0.75, 3.0, 0.5, 1.0

floats:
 	.float 0.0, 1.0

viewport:
	.float 0.002350, 0.004167 	@ dx dy

diffuse:
	.float 0.4, 0.7, 1.0
epsilon:
	.float 0.001

