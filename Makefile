NAME=intro

TOOLPREFIX=arm-linux-gnueabi

all: ${NAME}.bin
	qemu-system-arm -M vexpress-a9 -m 128M -kernel $^

%.elf: %.s
	${TOOLPREFIX}-as $^ -o $@
%.out: %.elf
	${TOOLPREFIX}-ld -Ttext=0x60010000 $^ -o $@
%.bin: %.out
	${TOOLPREFIX}-objcopy -O binary $^ $@
%.objdump: %.out
	${TOOLPREFIX}-objdump -d $^  > $@

clean:
	rm ${NAME}.bin

