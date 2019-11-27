#!/bin/bash
for j in $(seq 4 1 11)
do
	a=$((2 ** $j))
	make
	for i in $(seq 1 1 100)
	do
   		./check_memory
	done
	if [ $a -eq 16 ]
   	then
        	sed -i 's/#define N 16/#define N 32/g' "check_memory.c"
	fi
	if [ $a -eq 32 ]
   	then
        	sed -i 's/#define N 32/#define N 64/g' "check_memory.c"
	fi
        if [ $a -eq 64 ]
   	then
        	sed -i 's/#define N 64/#define N 128/g' "check_memory.c"
	fi
	if [ $a -eq 128 ]
   	then
        	sed -i 's/#define N 128/#define N 256/g' "check_memory.c"
	fi
	if [ $a -eq 256 ]
   	then
        	sed -i 's/#define N 256/#define N 512/g' "check_memory.c"
	fi
	if [ $a -eq 512 ]
   	then
        	sed -i 's/#define N 512/#define N 1024/g' "check_memory.c"
	fi
	if [ $a -eq 1024 ]
   	then
        	sed -i 's/#define N 1024/#define N 2048/g' "check_memory.c"
	fi

done

