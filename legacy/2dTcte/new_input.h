/*************** new_input.h  *************/
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define _GNU_SOURCE
#include <getopt.h>
#include <sys/types.h>



#define MAXINPUT 150
#define MAXM 10

#ifdef READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif


/* strings should have room for 80 characters */
/* or the length has to be stored somewhere? */

struct _input {
  char name[60];
  void *addr;
  char type;  /* 'i', 'd', 's', 'S'*/
  char tag[MAXM];  /* for menu and command line */
} linput[MAXINPUT];

int iinput = 0;


input_reset()
{
  iinput = 0;
}

input_insert(char * tag, char * name, void * addr, char type)
{
  strcpy(linput[iinput].name, name);
  linput[iinput].addr = addr;
  linput[iinput].type = type;
  strncpy(linput[iinput].tag, tag, MAXM);
  iinput ++;
}

input_redraw()
{
  int i;
  char s[100];
  
  printf("\n");
  for (i=0; i<iinput; i++) {
    printf("%*.*s (%s): ", MAXM, MAXM, linput[i].tag, linput[i].name);
    switch(linput[i].type){
    case 'i':
      printf("%d\n", * ((int *) linput[i].addr));
      break;
    case 'd':
      printf("%lg\n", *((double *) linput[i].addr));
      break;
    case 's':
      printf("'%s'\n", (char *) linput[i].addr);
      break;
    case 'S':
      input_interpolate(s,(char *) linput[i].addr);
      printf("'%s' ==> %s\n", (char *) linput[i].addr, s);
      break;
    }
  }
}


input_resume()
{
  int i;
  
  for (i=0; i<iinput; i++) {
    printf("%*.*s) %s: \t", MAXM,MAXM,linput[i].tag, linput[i].name);
    switch(linput[i].type){
    case 'i':
      printf("%d\n", * ((int *) linput[i].addr));
      break;
    case 'd':
      printf("%lg\n", *((double *) linput[i].addr));
      break;
    case 's':
    case 'S':
      printf("%s\n", (char *) linput[i].addr);
      break;
    }
  }
}

input_save(char * file) {
  int i;
  FILE *f;
  
  f = fopen(file, "w");
  print_header(f,"#");
  fprintf(f, "# parameters: \n");
  for (i=0; i<iinput; i++) {
    fprintf(f,"%*.*s : ", MAXM,MAXM,linput[i].tag);
    switch(linput[i].type){
    case 'i':
      fprintf(f,"%d  # %s\n", * ((int *) linput[i].addr), linput[i].name);
      break;
    case 'd':
      fprintf(f,"%lg # %s\n", *((double *) linput[i].addr), linput[i].name);
      break;
    case 's':
    case 'S':    
      fprintf(f,"%s # %s\n", (char *) linput[i].addr, linput[i].name);
      break;
    }
  }
  fclose(f);
}
char * trim(char * ss) {
	char * ws;
	/* eliminiamo gli spazi bianchi finali */
	ws=ss+strlen(ss)-1;
	while (isblank(*ws)){*(ws--)='\0';};
	ws=ss;
	/* eliminiamo gli spazi bianchi iniziali */
	while (isblank(*ws)){ws++;};
	return(ws);
}


int input_load(char * file) {
  int i;
  FILE *f;
  char ss[100];
	char * ws;
	char * ts;
  char s[100];
  char tag[100];
  int res;
  int line;
  
  f = fopen(file, "r");
  if (f==NULL) return(-1);
  for  (line=1; ; line++) {
    fgets(ss,80, f);
    if(feof(f)) break;
		/* eliminiamo i commenti */
		ws=strchr(ss, '#');
		if (ws) {
			*ws--='\0';
			while (isblank(*ws--)){};
		}
		if(strlen(ss)==0) continue;
		/*split su : */
		ws = strchr(ss, ':');
    if (!ws) {
      fprintf(stderr, "error in input file %s at line %d", file, line);
      exit(1);
    }  
		*ws++='\0';
		/* eliminiamo gli spazi bianchi */
		ws=trim(ws);
		ts=trim(ss);
    for (i=0; i<iinput; i++) { /* brutal linear search */
      if (strcmp(linput[i].tag, ts)==0) break;
    }
    
    if (i < iinput) {    
      switch(linput[i].type){
      case 'i':
        sscanf(ws,"%d", (int *) linput[i].addr);
        break;
      case 'd':
        sscanf(ws,"%lf", (double *) linput[i].addr);
        break;
      case 's':
      case 'S':
        sscanf(ws,"%s #", (char *) linput[i].addr);
        break;
      }
    } else {
      printf("unrecognized tag '%s' at line %d of file %s\n", ws,line,file);
    }
  }
  fclose(f);
}

char programname[100];

input_options(int argc, char ** argv, char * help)
{
  int interactive;
  int i,c;
  char * short_options= "?$";
  struct option long_options[MAXINPUT+1];
  int option_index;
  
	strcpy(programname, argv[0]);
  interactive = 0;
  for (i=0; i<iinput; i++) {
    long_options[i].name = linput[i].tag;
    long_options[i].has_arg = 1;
    long_options[i].flag = NULL;
    long_options[i].val = 0;
  }
  long_options[i].name = NULL;
  long_options[i].has_arg = 0;
  long_options[i].flag = NULL;
  long_options[i].val = 0;
  
  while ((c = getopt_long_only(argc, argv, short_options,
    long_options, &option_index)) != EOF) {
    if (c == '?') {
      printf("%s\n%s\ncompiled: %s on %s\n",SOURCE, help, DATE, HOSTNAME);
      printf("usage: %s [options]\n", argv[0]);
      printf("-? : this help\n");
      printf("-$ : interactive mode\n");
      for (i=0; i<iinput; i++) {
      	printf("--%0.*s <%s> (type %c)\n", MAXM,linput[i].tag, linput[i].name,
	      linput[i].type);
      }
			exit(0);
    } else if (c == '$') {
      interactive = 1;
    } else {
      i = option_index;
      if (i < iinput) {
	      switch(linput[i].type) {
	      case 'i':
	        *((int *) linput[i].addr) = atoi(optarg);
	        break;
	      case 'd':
	        *((double *) linput[i].addr) = atof(optarg);
	        break;
	      case 's':
        case 'S':
	        strcpy ((char *) linput[i].addr, optarg);
	        break;
	      }
      }
    }
  }
  return(interactive);
}



#ifdef READLINE
  static char * s = (char *)NULL;
  static char * s1 = (char *)NULL;
#endif


int input() {
  int i, c;
  int flag=1;
#ifdef READLINE
  char prompt[100];
#else
  char s[100];
  char s1[100];
#endif
  
  input_redraw();
  for (;;) {
#ifdef READLINE
   for (i=0; i<iinput; i++) { 
      sprintf(prompt, "%0.*s", MAXM, linput[i].tag);
      add_history(prompt);
    }
    strcpy(prompt, "\nchoose an option (. to accept) :> ");
    if (s1) {
      free(s1);
      s1 = (char *) NULL;
    } 
    s1=readline(prompt);
#else
    printf("\nchoose an option (. to accept) :> ");
    fgets(s1, 99, stdin);	
#endif
    if (s1[0] == '.') return 0;
    for (i=0; i<iinput; i++) { /* brutal linear search */
      if (strncmp(linput[i].tag, s1, MAXM)==0) break;
    }
    if (i < iinput) {
#ifdef READLINE
      switch(linput[i].type){
        case 'i':
          sprintf(prompt, "%d", * ((int *) linput[i].addr));
          break;
        case 'd':
          sprintf(prompt, "%lg", *((double *) linput[i].addr));
          break;
        case 's':
        case 'S':
          sprintf(prompt, "%s", (char *) linput[i].addr);
          break;
      }
      add_history(prompt);
      sprintf(prompt, "\n%s: ",linput[i].name);
      if (s) {
        free(s);
        s = (char *) NULL;
      }
      s=readline(prompt);
#else
      printf("\n%s: ",linput[i].name);
      fgets(s, 99, stdin);
#endif
      switch(linput[i].type){
        case 'i':
          sscanf(s,"%d", (int *) linput[i].addr);
          break;
        case 'd':
          sscanf(s,"%lf", (double *) linput[i].addr);
          break;
        case 's':                                                    
        case 'S':
          sscanf(s,"%s", (char *) linput[i].addr);
          break;
        default:
          break;
      }
      input_redraw(); 
    } else {
      printf("\nno such option\n");
      input_redraw();
    }                                                             
  }
}

/* input_interpolate prints formatted inputs into buf */
/* format is $5.2{var} or $5.2v if var name is 1 char */
/* isolated % characters are not allowed (use %% ) */

input_interpolate(char * buf, char * fmt) {
  char * tmp;
	char * i, * j, *z;
  char s[MAXM];
	int k;
  
  tmp = (char *) malloc(2*strlen(fmt)+100); /* mah...*/
  strcpy(buf, fmt);

  z=buf;
 	while (i=strchr(z, '$')) {
    z=i;
		for(; !(isalpha(*i) || *i=='{'); i++){};
    j=i;
    if (isalpha(*i)) {
      s[0] = *i;
      s[1]='\0';
    } else {
      i++;
  		for (k=0; *i!='}' && k<MAXM; k++,i++) {
        s[k]=*i;        
      }
      if(k<MAXM) {
        s[k] = '\0';
      }
    }
    for (k=0; k<iinput; k++) { /* brutal linear search */
      if (strncmp(linput[k].tag,s,MAXM)==0) break;
    }
    if (k>=iinput) { /* error, not found */
      z++;
      continue;  
    } else {
      *z='%';
      strcpy(j+1, i+1);
		  switch (linput[k].type) {
        case 'i' :
          *j = 'd';
          sprintf(tmp, buf,  *( (int*) linput[k].addr));
          break;
        case 'd' :
          *j = 'f';
          sprintf(tmp, buf,  *((double *) linput[k].addr));
          break;
        case 's' :
        case 'S' :
          *j = 's';
          sprintf(tmp, buf,  (char *) linput[k].addr);
          break;
      }
  		strcpy(buf, tmp);
    }
	}
  free(tmp);
  return(0);
}    

print_header(FILE *fout, char *comment) {
  char buf[128];
	time_t t;

	fprintf(fout, "%s program: %s\n",comment, programname);
	getcwd(buf, 128);
	fprintf(fout, "%s cwd: %s\n", comment, buf);
	gethostname(buf, 128);
	fprintf(fout,"%s hostname: %s\n", comment, buf);
	time(&t);
	fprintf(fout,"%s localtime: %s#\n", comment, ctime(&t));
}

print_inputs(FILE *fout, char *comment) {
  int i;
 
  print_header(fout, comment);  
	for (i=0; i<iinput; i++) {
    fprintf(fout,"%s %*.*s) %s: ", comment, MAXM,MAXM,linput[i].tag, 
      linput[i].name);
  	switch(linput[i].type){
    	 case 'i':
        	fprintf(fout, "%d\n",  *((int *)linput[i].addr)); 
        	break;
    	 case 'd':
        	fprintf(fout, "%lf\n", *((double*)linput[i].addr)); 
    	 break;
    	 case 's':
       case 'S':
        	fprintf(fout, "%s\n", (char *) linput[i].addr); 
    	 break;
    	 default:
    	 break;
  	}
	}
}

#if 0 
/* che fa? */
check_input(int i, char * str, char *comment)
{
 	char tstr[80];
  switch(linput[i].type){
     case 'i':
        sprintf(tstr, "%s %c) %s: %d\n", comment, linput[i].tag, 
           linput[i].name, *((int *)linput[i].addr)); 
        break;
     case 'd':
        sprintf(tstr, "%s %c) %s: %lf\n", comment, linput[i].tag, 
           linput[i].name, *((double*)linput[i].addr)); 
     break;
     case 's':
     case 'S':
        sprintf(tstr, "%s %c) %s: %s\n", comment, linput[i].tag, 
           linput[i].name, (char *) linput[i].addr); 
     break;
     default:
     break;
  }
	return(strcmp(str,tstr));
}
#endif

/* event manager (in the future...) */

int event_count=0;
int event_v = 0;

int input_event_setup(int v){
  event_count=0;
  event_v = v;
  if (v) {
    printf("%s", "0%       20%       40%       60%       80%      100%\n");
    printf("|----|----|----|----|----|----|----|----|----|----|\n ");
  }
  return 0;
}

int input_event(int n, int tot) {
  int i;
  int perc;
  
  perc = (n+1)*50/tot;
  if(perc > event_count) {
    if (event_v) {
      for(i=0; i<perc-event_count; i++) printf(".");
      if (perc == 50) {
        printf("\n");
      }
      fflush(stdout);
    }    
    event_count = perc;
  } 
  return 0;
}
  


