@include "/home/ehm/awkfiles/awklib.awk"
BEGIN{
  #
  # version 2.2 Fri Feb  6 16:34:33 CST 2015
  # -add switch for printing only atoms with atomic number Z
  #
  # version 2.1 Sun Feb  1 12:45:04 CST 2015
  # -add switch for not printing atoms with atomic number Z
  #
  # version 2.0, allow additional powder output types
  # default powderout=P (lattice vectors)
  # set powderout=L to use a,b,c, alph, beta, gamm
  # set molsym=1 to output brute force molecular symmetry finder.
  #

  # set on command line
  #
  #ignoreZ=0;
  #onlyZ=0;

  numZ=0;


  VERSION=2.2;
  CONVFMT="%.15g";
  r1[3];
  r2[3];
  r3[3];
  p[3000];
  debug=0;
  nat=0;
  count=0;

  if ( length(powderout)==0 ) powderout="P";
}
# note: this is coded to work with the output of
# icsd2poscar.awk, in which
# the first line should read: "Z: Z1 Z2 Z3 ..."
(NR==1){ for(i=2;i<NF+1;i++) { elem[i-1]=$i; ntyp=(i-1); } }

# if I change the first line of output back to
# Ca B H
# i.e., without the 'Z:' then this line works
# and looks cleaner and is easier to read
# (NR==1){ for(i=2;i<NF+1;i++) { elem[i]=$i; ntyp=i); } }

(NR==2){ scalefact=$1; }
(NR==3){ for(i=1;i<4;i++) { r1[i]=$i*scalefact; } }
(NR==4){ for(i=1;i<4;i++) { r2[i]=$i*scalefact; } }
(NR==5){ for(i=1;i<4;i++) { r3[i]=$i*scalefact; } }
(NR==6){ for(i=1;i<NF+1;i++) { numtyp[i]=$i; nat+=$i; } }
(NR>6 && NF>1) {
  p[3*count+0]=$1;
  p[3*count+1]=$2;
  p[3*count+2]=$3;
  count++;
}

END{
    
#  printf("numtypes = %d\n",ntyp);
# for(j=1;j<ntyp+1;j++) { printf("%d ",numtyp[j]); }
    if ( ignoreZ > 0 ){
	for(j=1;j<ntyp+1;j++){if (elem[j]==ignoreZ) {numZ+=numtyp[j]} }
	printf("%d\n",nat-numZ);
    } else if ( onlyZ > 0 ) {
	for(j=1;j<ntyp+1;j++){if (elem[j]==onlyZ) {numZ+=numtyp[j]} }
	printf("%d\n",numZ);
    } else {
	printf("%d\n",nat);
    }
    if ( molsym != 1 && powderout=="P" ){
	printf("%20.15f%20.15f%20.15f\n",r1[1],r1[2],r1[3]);
	printf("%20.15f%20.15f%20.15f\n",r2[1],r2[2],r2[3]);
	printf("%20.15f%20.15f%20.15f\n",r3[1],r3[2],r3[3]);
	printf("\n");
    }
    if ( powderout=="L" ){
	R1[0]=r1[1];      R1[1]=r1[2];      R1[2]=r1[3];
	R2[0]=r2[1];      R2[1]=r2[2];      R2[2]=r2[3];
	R3[0]=r3[1];      R3[1]=r3[2];      R3[2]=r3[3];
	latvec2abc(R1,R2,R3,"L");
    }
    count=0;
    for(j=1;j<ntyp+1;j++) {
	elm=elem[j];
	if ( elm == ignoreZ ) { count=count+numtyp[j]; continue; }
	if ( onlyZ && elm != onlyZ ) { count=count+numtyp[j]; continue; }
	for(i=0;i<numtyp[j];i++){
	    at=count;
	    if ( molsym==1 ){
		printf("%5s%20.15f%20.15f%20.15f\n",elm,p[3*at+0],p[3*at+1],p[3*at+2]);
	    } else {
		printf("%20.15f%20.15f%20.15f%5s\n",p[3*at+0],p[3*at+1],p[3*at+2],elm);
	    }
	    count++;
	}
	printf("\n");
    }
    
}
