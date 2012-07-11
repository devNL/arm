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
	call max
	vmov s2, s1
	call min
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

# 
normalize:

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
	vmov s5, s0 #tmp=bla
	vsub.f32 s5, s3 #tmp-=tx
	vmov s0, s5 #arg0=bla-tx
	vmov s1, s2 #arg1=pz
	bl length2
	vsub.f32 s0, s4
bx lr

# args: px, py, pz
dist:

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

normal_at_point:

bx lr


torus:
	.float 4.0, 0.0, 9.0

box:
	.float -3.0, 0.0, 10.0
