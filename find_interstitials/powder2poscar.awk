# Calls to awklib.awk
@include "/home/ehm/awkfiles/awklib.awk"
BEGIN{

  CONVFMT="%.15g";

  # set on command line
  # PFT = 1,2,3  powder file type (see powder -F for list)

  # debug=0;

  r1[3];
  r2[3];
  r3[3];
  p[3000];
  Z[1000];

  nat=0;
  count=0;
  record=0;

  # allow 1000 types of different atoms
  n_type=0;
  type[1000];
  ntype[1000];
  for(i=0;i<1000;i++){ type[i]= -1; ntype[i]=0; }

}

# note: this is coded to work with the output of
# icsd2poscar.awk, in which
# the first line should read: "Z: Z1 Z2 Z3 ..."

# the only line of the powder file having one field
# should be the number of atoms

##########################
#  Functions
##########################



##########################
#  pattern action rules
##########################
(NF==1){num_atoms=$1;}

(NF==3 && record==0 && PFT==3){ for(i=1;i<4;i++){ r1[i]=$i; };record=NR;}
(NF==3 && NR==record+1 && PFT==3){ for(i=1;i<4;i++) { r2[i]=$i; }; }
(NF==3 && NR==record+2 && PFT==3){ for(i=1;i<4;i++) { r3[i]=$i; }; }

(NF==3 && record==0 && PFT==2){alat=$1;blat=$2;clat=$3;record=NR;}
(NF==3 && NR==record+1 && PFT==2){alph=$1;beta=$2;gamm=$3;}

(NF==4){
  p[3*count+0]=$1;
  p[3*count+1]=$2;
  p[3*count+2]=$3;
  Z[count]=$4;
  count++;
}


END{

  if ( ! PFT ){
    printf("Must set PFT, powder file type. 1,2, or 3\n");
    printf("See powder -F for file types.\n");
    exit;
  }
# note: this is coded to work with the output of
# icsd2poscar.awk, in which
# the first line should read: "Z: Z1 Z2 Z3 ..."

  for(i=0;i<count;i++){
    if ( debug) printf("Z[%d] = %d\n",i,Z[i]);

    match_flag=0;
    for(j=0;j<n_type;j++){
      if ( debug ) printf("j=%d\n",j);
      if ( type[j] == Z[i] ) {
	if ( debug ) printf("found match to type %d\n",j);
	ntype[j]++;
	match_flag=1;
	break;
      }
    }
    if ( match_flag==0 ) {
      if ( debug) printf("setting type[%d] = %d\n",n_type,Z[i]);
      type[n_type]=Z[i];
      n_type++;
    }

  }

  if ( debug ) {
    printf("Found n_type = %d\n",n_type);
    for(i=0;i<n_type;i++){
      printf(" type[%d] = %2d\n",i,type[i]);
      printf("ntype[%d] = %2d\n",i,ntype[i]+1);
    }
  }

  if ( debug > 0 ){
    printf("alat blat clat = %f  %f  %f\n",alat,blat,clat);
    printf("alph beta gamm = %f  %f  %f\n",alph,beta,gamm);
  }

  printf("Z: ");
  for(i=0;i<n_type;i++){ printf("%d ",type[i]); }
  printf("\n");
  if ( PFT==2 ){ print_lats(alat,blat,clat,alph,beta,gamm);}
  if ( PFT==3 ){
    printf("1.0 lattice constant\n");
    printf("%20.15f%20.15f%20.15f\n",r1[1],r1[2],r1[3]);
    printf("%20.15f%20.15f%20.15f\n",r2[1],r2[2],r2[3]);
    printf("%20.15f%20.15f%20.15f\n",r3[1],r3[2],r3[3]);
  }
  for(i=0;i<n_type;i++){ printf("%5d ",ntype[i]+1); }
  printf("\n");
  printf("Direct\n");

  for(i=0;i<n_type;i++){
    for(j=0;j<count;j++){
      if ( Z[j] == type[i] ){
	printf("%20.15f%20.15f%20.15f\n",p[3*j+0],p[3*j+1],p[3*j+2]);
      }
    }
    printf("\n");
  }
  
}
