%.e: %.c arrays.h sol.h  new_input.h tools.h 
	gcc $< -o $@  -DREADLINE -DREUSE -DSOURCE=\"$<\" -DDATE=\""`date`"\" -DHOSTNAME=\"`hostname`\"  -lreadline -lncurses -lm    -lpq -I /usr/include/postgresql/ 

#-lpq -I /usr/include/postgresql/
