.arch armv7-a
.fpu neon
.syntax unified

.global _start
_start:

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



mov r0, $0
mov r3, $0x12c000

.redbars:
mov r1, r0
and r1, r1, $0xFF
str r1, [r2, r0]
add r0, r0, $4
cmp r0, r3
bne .redbars

hcf: B hcf

