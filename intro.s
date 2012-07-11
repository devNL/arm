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

# args x,a,b
# return (x<a)?a:(x>b)?b:x;
clamp:
bx lr

# args: x,a
# return (x<a)?x:a;
min:
	# if (r0 < r1) return r0
	cmp r0, r1
	blt end_min
	# else return r1
	mov r0, r1 
end_min:
	bx lr

# args: x,a
# return (x>a)?x:a;
max:
bx lr

# args: x,y
# return sqrtf(x*x+y*y);
length2:
bx lr

# args: x,y,z
# return sqrtf(x*x+y*y+z*z);
length3:
bx lr

# 
normalize:
bx lr

mix:
bx lr

udroundbox:
bx lr

sdtorus:
bx lr

dist:
bx lr

dot:
bx lr

normal_at_point:
bx lr
