all:  average.e

%.e: %.c arrays.h  new_input.h tools.h 
	gcc $< -o $@  -DREADLINE -DREUSE -DSOURCE=\"$<\" -DDATE=\""`date`"\" -DHOSTNAME=\"`hostname`\"  -lreadline -lncurses -lm -O3 -p -g  
