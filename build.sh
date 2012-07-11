#!/bin/sh

arm-linux-as -o intro.elf intro.s
arm-linux-objcopy -O binary intro.elf intro.bin
