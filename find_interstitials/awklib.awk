BEGIN{

  # to use add: @include "/home/ehm/awkfiles/awklib.awk"
  # to the first line of your awk file, and call with igawk
  # OR set your AWKPATH to the directory containing this file

  version="1.9.5"

  # v 1.9.5 Mon Jun 26 07:21:00 CDT 2017
  #         add function tag_value
  #
  # v 1.9.4 Thu June 08 9:02 CST 2017
  #         move some direct to cart code into awklib
  #         for use in another script to convert XDATCAR
  #         info to an xbsa movie file
  #
  # v 1.9.3 Fri Feb 13 16:32:17 CST 2015
  #         add 3x3 M on vector multiply
  #
  # v 1.9.2 Tue 19 Nov 2013
  #         add to latvec2abc a powderout=L option
  #
  # v 1.9.1 Wed May  1 15:59:58 CDT 2013
  #         fix iabs and fabs functions
  #
  # v 1.9 Mon Nov  1 14:29:03 CDT 2010
  #       add iabs and fabs functions
  #
  # v 1.8 Sun Feb 28 22:11:11 CST 2010
  #       add new function wttoz for apos.dat to POSCAR
  #
  # v 1.7 Tue Jun 30 10:01:23 CDT 2009
  #       add new function for checking whether a cell
  #       is cubic or not.
  #
  # v 1.6 Wed Jul 30 13:39:08 CDT 2008
  #       fix some of the vadd, vsub, codes
  #       It seems that vmag() works and that gawk can pass
  #       an array to a function and use it with a different
  #       name *inside* the function, and return a scalar.
  #       But if one tries to pass an array name back with
  #       return, it cannot do it.  One has
  #       to alter the array inside the function.  So vadd and vsub
  #       are probably useless, the calling script has to use v3.
  #
  # v 1.5 Wed Jan 16 16:19:22 CST 2008
  #       add output specific to powder -L option
  #
  # v 1.4 Tue Jul 31 15:41:50 CDT 2007
  #       added atomic weights for density in contcar_info
  #
  # v 1.3 Fri Jun 22 13:48:26 PDT 2007
  #       added some vector ops
  #
  # v 1.2 Wed Jun 20 17:19:16 PDT 2007
  #       added function ztoelm from poscar2vis.awk
  #
  # v 1.1 Tue Apr  3 09:12:47 PDT 2007
  #       fixed triple product bug


  CONVFMT="%.15g";
  PI =  4.0*atan2(1.0,1.0);
  D2R = PI/180.0;
  R2D = 180.0/PI;
}

function iabs (n) {
  if ( n < 0 ) return (-n);
  else return (n);
}

function fabs (a) {
  if ( a < 0.0 ) return (-a);
  else return (a);
}

function tol (a,b,t) {
  if ( ( (a+t)>b ) && ( (a-t)<b ) ) return 1;
  else return 0;
}

function vmag (a) {
  return sqrt( a[0]*a[0]+a[1]*a[1]+a[2]*a[2] );
}

function vnorm (a) {
  amag = vmag(a);
  for(indext1=0;indext1<3;indext1++) norm[indext1] = a[indext1]/amag;
  return;
}

function vadd (v1,v2) {
  v3[0]=v1[0]+v2[0];
  v3[1]=v1[1]+v2[1];
  v3[2]=v1[2]+v2[2];
  return;
}

function vsub (v1,v2) {
  v3[0]=v1[0]-v2[0];
  v3[1]=v1[1]-v2[1];
  v3[2]=v1[2]-v2[2];
  return;
}

function vdotprod (v1,v2) {
  return ( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] );
}

function dotprod(v1,v2,v3,w1,w2,w3){
  return v1*w1+v2*w2+v3*w3;
}

function angle(bb,cc){
  aa = sqrt(cc*cc - bb*bb);
  theta = R2D * atan2(aa,bb);
  return theta;
}

function print_matrix_3x3(M){
  
  printf("Matrix:\n");
  printf("%20.15f%20.15f%20.15f\n",M[0],M[1],M[2]);
  printf("%20.15f%20.15f%20.15f\n",M[3],M[4],M[5]);
  printf("%20.15f%20.15f%20.15f\n",M[6],M[7],M[8]);

  return;
}

function matrix_33_vec_mult(M,vin,vout){
    vout[0] = M[0]*vin[0] + M[1]*vin[1] + M[2]*vin[2];
    vout[1] = M[3]*vin[0] + M[4]*vin[1] + M[5]*vin[2];
    vout[2] = M[6]*vin[0] + M[7]*vin[1] + M[8]*vin[2];
}

#
# a[0] a[1] a[2] a[3] a[4] a[5] a[6] a[7] a[8]
# a11  a12  a13  a21  a22  a23  a31  a32  a33
#
function det(M){
  return ( M[0]*(M[8]*M[4]-M[7]*M[5]) - M[3]*(M[8]*M[1]-M[7]*M[2]) + M[6]*(M[5]*M[1]-M[4]*M[2]) );
}

function transpose(M,Mt){

  Mt[0] = M[0];
  Mt[1] = M[3];
  Mt[2] = M[6];
  Mt[3] = M[1];
  Mt[4] = M[4];
  Mt[5] = M[7];
  Mt[6] = M[2];
  Mt[7] = M[5];
  Mt[8] = M[8];

  return;
}

function inverse3x3(M,Minv){
  
  D = det(M);
  if ( tol(D,0.0,1e-15 ) ) {
    printf("Determinant is zero!! --> exiting.\n");
    exit;
  }
  Minv[0] =  (M[8]*M[4]-M[7]*M[5]) / D;
  Minv[1] = -(M[8]*M[1]-M[7]*M[2]) / D;
  Minv[2] =  (M[5]*M[1]-M[4]*M[2]) / D;
  Minv[3] = -(M[8]*M[3]-M[6]*M[5]) / D;
  Minv[4] =  (M[8]*M[0]-M[6]*M[2]) / D;
  Minv[5] = -(M[5]*M[0]-M[3]*M[2]) / D;
  Minv[6] =  (M[7]*M[3]-M[6]*M[4]) / D;
  Minv[7] = -(M[7]*M[0]-M[6]*M[1]) / D;
  Minv[8] =  (M[4]*M[0]-M[3]*M[1]) / D;
  
  return;
}

function iscubic(M){
  CUBTOL=1e-5;
  if ( tol(M[0],M[4],CUBTOL) && tol(M[0],M[8],CUBTOL) &&
       tol(M[1],0.0,CUBTOL) && tol(M[2],0.0,CUBTOL) &&
       tol(M[3],0.0,CUBTOL) && tol(M[5],0.0,CUBTOL) &&
       tol(M[6],0.0,CUBTOL) && tol(M[7],0.0,CUBTOL) ){
    return(M[0]);
  } else return(0);
}

function latvec2abc(a,b,c,powderout){

  A = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
  B = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
  C = sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] );

  triple = c[0]*(a[1]*b[2]-a[2]*b[1]) + c[1]*(a[2]*b[0]-a[0]*b[2]) + c[2]*(a[0]*b[1]-a[1]*b[0]);

  # there is no acos in awk, have to find with atan2
  gm_b = dotprod(a[0],a[1],a[2],b[0],b[1],b[2]);
  gm_c = (A*B);
  if ( debug>0 ){printf("gm_b = %f\n",gm_b);}
  gamma = angle(gm_b,gm_c);

  al_b = dotprod(c[0],c[1],c[2],b[0],b[1],b[2]);
  al_c = (B*C);
  alpha = angle(al_b,al_c);

  bt_b = dotprod(c[0],c[1],c[2],a[0],a[1],a[2]);
  bt_c = (C*A);
  beta = angle(bt_b,bt_c);
  if ( powderout==0 ){
      printf("a = %15.10f%15.10f%15.10f\n",a[0],a[1],a[2]);
      printf("b = %15.10f%15.10f%15.10f\n",b[0],b[1],b[2]);
      printf("c = %15.10f%15.10f%15.10f\n",c[0],c[1],c[2]);
      printf("\n");
      printf("A = %15.10f\n",A);
      printf("B = %15.10f\n",B);
      printf("C = %15.10f\n",C);
      if ( triple < 0.0 ) {
	  printf("!!! WARN !!!: (A x B . C) is negative. Check cell vectors!\n");
	  triple *= -1;
      }
      printf("Vol = %15.10f\n",triple);
      printf("Vol = %15.10f  (if it were an orthorhombic cell)\n",A*B*C);
      
      printf("\n");
      printf("alpha = %15.10f\n",alpha);
      printf("beta  = %15.10f\n",beta);
      printf("gamma = %15.10f\n",gamma);
      printf("\n");
      printf("orth  = %.5f\n",triple/(A*B*C));
      
      printf("Output for pack restart file header:\n");
      printf("*r %20.10f%20.10f%20.10f\n",A,B,C);
      printf("*r %20.10f%20.10f%20.10f\n",alpha*D2R,beta*D2R,gamma*D2R);
      
      printf("\nOutput for powder -L option conventional cell\n");
      printf("*p %15.10f%15.10f%15.10f\n",A,B,C);
      printf("*p %15.10f%15.10f%15.10f\n\n",alpha,beta,gamma);
  }
  if ( powderout=="L" ){
      printf("%15.10f%15.10f%15.10f\n",A,B,C);
      printf("%15.10f%15.10f%15.10f\n",alpha,beta,gamma);      
  }
  
}

function cellvolume(a,b,c){

  A = sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
  B = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
  C = sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] );

  triple = c[0]*(a[1]*b[2]-a[2]*b[1]) + c[1]*(a[2]*b[0]-a[0]*b[2]) + c[2]*(a[0]*b[1]-a[1]*b[0]);
  if ( triple < 0.0 ) { triple *= -1; }

#  printf("a = %15.10f%15.10f%15.10f\n",a[0],a[1],a[2]);
#  printf("b = %15.10f%15.10f%15.10f\n",b[0],b[1],b[2]);
#  printf("c = %15.10f%15.10f%15.10f\n",c[0],c[1],c[2]);
#  printf("\n");
#  printf("A = %15.10f\n",A);
#  printf("B = %15.10f\n",B);
#  printf("C = %15.10f\n",C);
#  printf("Vol = %15.10f\n",triple);

  return( triple );
}


function get_Z( ELMT ) {
  if ( ELMT=="H" )  return   1;
  if ( ELMT=="He" ) return   2;
  if ( ELMT=="Li" ) return   3;
  if ( ELMT=="Be" ) return   4;
  if ( ELMT=="B" )  return   5;
  if ( ELMT=="C" )  return   6;
  if ( ELMT=="N" )  return   7;
  if ( ELMT=="O" )  return   8;
  if ( ELMT=="F" )  return   9;
  if ( ELMT=="Ne" ) return  10;
  if ( ELMT=="Na" ) return  11;
  if ( ELMT=="Mg" ) return  12;
  if ( ELMT=="Al" ) return  13;
  if ( ELMT=="Si" ) return  14;
  if ( ELMT=="P" )  return  15;
  if ( ELMT=="S" )  return  16;
  if ( ELMT=="Cl" ) return  17;
  if ( ELMT=="Ar" ) return  18;
  if ( ELMT=="K" )  return  19;
  if ( ELMT=="Ca" ) return  20;
  if ( ELMT=="Sc" ) return  21;
  if ( ELMT=="Ti" ) return  22;
  if ( ELMT=="V" )  return  23;
  if ( ELMT=="Cr" ) return  24;
  if ( ELMT=="Mn" ) return  25;
  if ( ELMT=="Fe" ) return  26;
  if ( ELMT=="Co" ) return  27;
  if ( ELMT=="Ni" ) return  28;
  if ( ELMT=="Cu" ) return  29;
  if ( ELMT=="Zn" ) return  30;
  if ( ELMT=="Ga" ) return  31;
  if ( ELMT=="Ge" ) return  32;
  if ( ELMT=="As" ) return  33;
  if ( ELMT=="Se" ) return  34;
  if ( ELMT=="Br" ) return  35;
  if ( ELMT=="Kr" ) return  36;
  if ( ELMT=="Rb" ) return  37;
  if ( ELMT=="Sr" ) return  38;
  if ( ELMT=="Y" )  return  39;
  if ( ELMT=="Zr" ) return  40;
  if ( ELMT=="Nb" ) return  41;
  if ( ELMT=="Mo" ) return  42;
  if ( ELMT=="Tc" ) return  43;
  if ( ELMT=="Ru" ) return  44;
  if ( ELMT=="Rh" ) return  45;
  if ( ELMT=="Pd" ) return  46;
  if ( ELMT=="Ag" ) return  47;
  if ( ELMT=="Cd" ) return  48;
  if ( ELMT=="In" ) return  49;
  if ( ELMT=="Sn" ) return  50;
  if ( ELMT=="Sb" ) return  51;
  if ( ELMT=="Te" ) return  52;
  if ( ELMT=="I" )  return  53;
  if ( ELMT=="Xe" ) return  54;
  if ( ELMT=="Cs" ) return  55;
  if ( ELMT=="Ba" ) return  56;
  if ( ELMT=="La" ) return  57;
  if ( ELMT=="Ce" ) return  58;
  if ( ELMT=="Pr" ) return  59;
  if ( ELMT=="Nd" ) return  60;
  if ( ELMT=="Pm" ) return  61;
  if ( ELMT=="Sm" ) return  62;
  if ( ELMT=="Eu" ) return  63;
  if ( ELMT=="Gd" ) return  64;
  if ( ELMT=="Tb" ) return  65;
  if ( ELMT=="Dy" ) return  66;
  if ( ELMT=="Ho" ) return  67;
  if ( ELMT=="Er" ) return  68;
  if ( ELMT=="Tm" ) return  69;
  if ( ELMT=="Yb" ) return  70;
  if ( ELMT=="Lu" ) return  71;
  if ( ELMT=="Hf" ) return  72;
  if ( ELMT=="Ta" ) return  73;
  if ( ELMT=="W" )  return  74;
  if ( ELMT=="Re" ) return  75;
  if ( ELMT=="Os" ) return  76;
  if ( ELMT=="Ir" ) return  77;
  if ( ELMT=="Pt" ) return  78;
  if ( ELMT=="Au" ) return  79;
  if ( ELMT=="Hg" ) return  80;
  if ( ELMT=="Tl" ) return  81;
  if ( ELMT=="Pb" ) return  82;
  if ( ELMT=="Bi" ) return  83;
  if ( ELMT=="Po" ) return  84;
  if ( ELMT=="At" ) return  85;
  if ( ELMT=="Rn" ) return  86;
  if ( ELMT=="Fr" ) return  87;
  if ( ELMT=="Ra" ) return  88;
  if ( ELMT=="Ac" ) return  89;
  if ( ELMT=="Th" ) return  90;
  if ( ELMT=="Pa" ) return  91;
  if ( ELMT=="U" )  return  92;
  if ( ELMT=="Np" ) return  93;
  if ( ELMT=="Pu" ) return  94;
  if ( ELMT=="Am" ) return  95;
  if ( ELMT=="Cm" ) return  96;
  if ( ELMT=="Bk" ) return  97;
  if ( ELMT=="Cf" ) return  98;
  if ( ELMT=="Es" ) return  99;
  if ( ELMT=="Fm" ) return 100;
  if ( ELMT=="Md" ) return 101;
  if ( ELMT=="No" ) return 102;
  if ( ELMT=="Lr" ) return 103;
  if ( ELMT=="Rf" ) return 104;
  if ( ELMT=="Db" ) return 105;
  if ( ELMT=="Sg" ) return 106;
  if ( ELMT=="Bh" ) return 107;
  if ( ELMT=="Hs" ) return 108;
  if ( ELMT=="Mt" ) return 109;
}

function ztoelm( z ) {
  if ( z==1   ) return( "H"  );
  if ( z==2   ) return( "He" );
  if ( z==3   ) return( "Li" );
  if ( z==4   ) return( "Be" );
  if ( z==5   ) return( "B"  );
  if ( z==6   ) return( "C"  );
  if ( z==7   ) return( "N"  );
  if ( z==8   ) return( "O"  );
  if ( z==9   ) return( "F"  );
  if ( z==10  ) return( "Ne" );
  if ( z==11  ) return( "Na" );
  if ( z==12  ) return( "Mg" );
  if ( z==13  ) return( "Al" );
  if ( z==14  ) return( "Si" );
  if ( z==15  ) return( "P"  );
  if ( z==16  ) return( "S"  );
  if ( z==17  ) return( "Cl" );
  if ( z==18  ) return( "Ar" );
  if ( z==19  ) return( "K"  );
  if ( z==20  ) return( "Ca" );
  if ( z==21  ) return( "Sc" );
  if ( z==22  ) return( "Ti" );
  if ( z==23  ) return( "V"  );
  if ( z==24  ) return( "Cr" );
  if ( z==25  ) return( "Mn" );
  if ( z==26  ) return( "Fe" );
  if ( z==27  ) return( "Co" );
  if ( z==28  ) return( "Ni" );
  if ( z==29  ) return( "Cu" );
  if ( z==30  ) return( "Zn" );
  if ( z==31  ) return( "Ga" );
  if ( z==32  ) return( "Ge" );
  if ( z==33  ) return( "As" );
  if ( z==34  ) return( "Se" );
  if ( z==35  ) return( "Br" );
  if ( z==36  ) return( "Kr" );
  if ( z==37  ) return( "Rb" );
  if ( z==38  ) return( "Sr" );
  if ( z==39  ) return( "Y"  );
  if ( z==40  ) return( "Zr" );
  if ( z==41  ) return( "Nb" );
  if ( z==42  ) return( "Mo" );
  if ( z==43  ) return( "Tc" );
  if ( z==44  ) return( "Ru" );
  if ( z==45  ) return( "Rh" );
  if ( z==46  ) return( "Pd" );
  if ( z==47  ) return( "Ag" );
  if ( z==48  ) return( "Cd" );
  if ( z==49  ) return( "In" );
  if ( z==50  ) return( "Sn" );
  if ( z==51  ) return( "Sb" );
  if ( z==52  ) return( "Te" );
  if ( z==53  ) return( "I"  );
  if ( z==54  ) return( "Xe" );
  if ( z==55  ) return( "Cs" );
  if ( z==56  ) return( "Ba" );
  if ( z==57  ) return( "La" );
  if ( z==58  ) return( "Ce" );
  if ( z==59  ) return( "Pr" );
  if ( z==60  ) return( "Nd" );
  if ( z==61  ) return( "Pm" );
  if ( z==62  ) return( "Sm" );
  if ( z==63  ) return( "Eu" );
  if ( z==64  ) return( "Gd" );
  if ( z==65  ) return( "Tb" );
  if ( z==66  ) return( "Dy" );
  if ( z==67  ) return( "Ho" );
  if ( z==68  ) return( "Er" );
  if ( z==69  ) return( "Tm" );
  if ( z==70  ) return( "Yb" );
  if ( z==71  ) return( "Lu" );
  if ( z==72  ) return( "Hf" );
  if ( z==73  ) return( "Ta" );
  if ( z==74  ) return( "W"  );
  if ( z==75  ) return( "Re" );
  if ( z==76  ) return( "Os" );
  if ( z==77  ) return( "Ir" );
  if ( z==78  ) return( "Pt" );
  if ( z==79  ) return( "Au" );
  if ( z==80  ) return( "Hg" );
  if ( z==81  ) return( "Tl" );
  if ( z==82  ) return( "Pb" );
  if ( z==83  ) return( "Bi" );
  if ( z==84  ) return( "Po" );
  if ( z==85  ) return( "At" );
  if ( z==86  ) return( "Rn" );
  if ( z==87  ) return( "Fr" );
  if ( z==88  ) return( "Ra" );
  if ( z==89  ) return( "Ac" );
  if ( z==90  ) return( "Th" );
  if ( z==91  ) return( "Pa" );
  if ( z==92  ) return( "U " );
  if ( z==93  ) return( "Np" );
  if ( z==94  ) return( "Pu" );
  if ( z==95  ) return( "Am" );
  if ( z==96  ) return( "Cm" );
  if ( z==97  ) return( "Bk" );
  if ( z==98  ) return( "Cf" );
  if ( z==99  ) return( "Es" );
  if ( z==100 ) return( "Fm" );
  if ( z==101 ) return( "Md" );
  if ( z==102 ) return( "No" );
  if ( z==103 ) return( "Lr" );
}

function ztowt( z ) {
  if ( z==1   ) return( 1.01  );
  if ( z==2   ) return( 4.00 );
  if ( z==3   ) return( 6.94 );
  if ( z==4   ) return( 9.01 );
  if ( z==5   ) return( 10.81  );
  if ( z==6   ) return( 12.01  );
  if ( z==7   ) return( 14.01  );
  if ( z==8   ) return( 16.00  );
  if ( z==9   ) return( 19.00  );
  if ( z==10  ) return( 20.18 );
  if ( z==11  ) return( 22.99 );
  if ( z==12  ) return( 24.30 );
  if ( z==13  ) return( 26.98 );
  if ( z==14  ) return( 28.09 );
  if ( z==15  ) return( 30.97  );
  if ( z==16  ) return( 32.06  );
  if ( z==17  ) return( 35.45 );
  if ( z==18  ) return( 39.95 );
  if ( z==19  ) return( 39.09  );
  if ( z==20  ) return( 40.08 );
  if ( z==21  ) return( 44.96 );
  if ( z==22  ) return( 47.90 );
  if ( z==23  ) return( 50.94  );
  if ( z==24  ) return( 52.00 );
  if ( z==25  ) return( 54.94 );
  if ( z==26  ) return( 55.85 );
  if ( z==27  ) return( 58.93 );
  if ( z==28  ) return( 58.71 );
  if ( z==29  ) return( 63.55 );
  if ( z==30  ) return( 65.38 );
  if ( z==31  ) return( 69.72 );
  if ( z==32  ) return( 72.59 );
  if ( z==33  ) return( 74.92 );
  if ( z==34  ) return( 78.96 );
  if ( z==35  ) return( 79.91 );
  if ( z==36  ) return( 83.80 );
  if ( z==37  ) return( 85.47 );
  if ( z==38  ) return( 87.62 );
  if ( z==39  ) return( 88.91  );
  if ( z==40  ) return( 91.22 );
  if ( z==41  ) return( 92.91 );
  if ( z==42  ) return( 95.94 );
  if ( z==43  ) return( 98.91 );
  if ( z==44  ) return( 101.07 );
  if ( z==45  ) return( 102.90 );
  if ( z==46  ) return( 106.40 );
  if ( z==47  ) return( 107.87 );
  if ( z==48  ) return( 112.40 );
  if ( z==49  ) return( 114.82 );
  if ( z==50  ) return( 118.69 );
  if ( z==51  ) return( 121.75 );
  if ( z==52  ) return( 127.60 );
  if ( z==53  ) return( 126.90  );
  if ( z==54  ) return( 131.30 );
  if ( z==55  ) return( 132.91 );
  if ( z==56  ) return( 137.34 );
  if ( z==57  ) return( 138.91 );
  if ( z==58  ) return( 140.12 );
  if ( z==59  ) return( 140.91 );
  if ( z==60  ) return( 144.24 );
  if ( z==61  ) return( 145.00 );
  if ( z==62  ) return( 150.35 );
  if ( z==63  ) return( 151.96 );
  if ( z==64  ) return( 157.25 );
  if ( z==65  ) return( 158.92 );
  if ( z==66  ) return( 162.50 );
  if ( z==67  ) return( 164.93 );
  if ( z==68  ) return( 167.26 );
  if ( z==69  ) return( 168.93 );
  if ( z==70  ) return( 173.04 );
  if ( z==71  ) return( 174.97 );
  if ( z==72  ) return( 178.49 );
  if ( z==73  ) return( 180.95 );
  if ( z==74  ) return( 183.85  );
  if ( z==75  ) return( 186.20 );
  if ( z==76  ) return( 190.20 );
  if ( z==77  ) return( 192.22 );
  if ( z==78  ) return( 195.09 );
  if ( z==79  ) return( 196.97 );
  if ( z==80  ) return( 200.59 );
  if ( z==81  ) return( 204.37 );
  if ( z==82  ) return( 207.20 );
  if ( z==83  ) return( 208.98 );
  if ( z==84  ) return( 210.00 );
  if ( z==85  ) return( 210.00 );
  if ( z==86  ) return( 222.00 );
  if ( z==87  ) return( 223.00 );
  if ( z==88  ) return( 226.02 );
  if ( z==89  ) return( 227.00 );
  if ( z==90  ) return( 232.04 );
  if ( z==91  ) return( 231.00 );
  if ( z==92  ) return( 238.03 );
  if ( z==93  ) return( 237.00 );
  if ( z==94  ) return( 244.00 );
  if ( z==95  ) return( 243.00 );
  if ( z==96  ) return( "Cm" );
  if ( z==97  ) return( "Bk" );
  if ( z==98  ) return( "Cf" );
  if ( z==99  ) return( "Es" );
  if ( z==100 ) return( "Fm" );
  if ( z==101 ) return( "Md" );
  if ( z==102 ) return( "No" );
  if ( z==103 ) return( "Lr" );
}


function wttoz( wt ) {
  if ( tol(wt,1.01,0.1)==1 ) return( 1 );
  if ( tol(wt,4.00,0.1)==1 ) return( 2 );
  if ( tol(wt,6.94,0.1)==1 ) return( 3 );
  if ( tol(wt,9.01,0.1)==1 ) return( 4 );
  if ( tol(wt,10.81,0.1)==1 ) return( 5 );
  if ( tol(wt,12.01,0.1)==1 ) return( 6 );
  if ( tol(wt,14.01,0.1)==1 ) return( 7 );
  if ( tol(wt,16.00,0.1)==1 ) return( 8 );
  if ( tol(wt,19.00,0.1)==1 ) return( 9 );
  if ( tol(wt,20.18,0.1)==1 ) return( 10 );
  if ( tol(wt,22.99,0.1)==1 ) return( 11 );
  if ( tol(wt,24.30,0.1)==1 ) return( 12 );
  if ( tol(wt,26.98,0.1)==1 ) return( 13 );
  if ( tol(wt,28.09,0.1)==1 ) return( 14 );
  if ( tol(wt,30.97,0.1)==1 ) return( 15 );
  if ( tol(wt,32.06,0.1)==1 ) return( 16 );
  if ( tol(wt,35.45,0.1)==1 ) return( 17 );
  if ( tol(wt,39.95,0.1)==1 ) return( 18 );
  if ( tol(wt,39.09,0.1)==1 ) return( 19 );
  if ( tol(wt,40.08,0.1)==1 ) return( 20 );
  if ( tol(wt,44.96,0.1)==1 ) return( 21 );
  if ( tol(wt,47.90,0.1)==1 ) return( 22 );
  if ( tol(wt,50.94,0.1)==1 ) return( 23 );
  if ( tol(wt,52.00,0.1)==1 ) return( 24 );
  if ( tol(wt,54.94,0.1)==1 ) return( 25 );
  if ( tol(wt,55.85,0.1)==1 ) return( 26 );
  if ( tol(wt,58.93,0.1)==1 ) return( 27 );
  if ( tol(wt,58.71,0.1)==1 ) return( 28 );
  if ( tol(wt,63.55,0.1)==1 ) return( 29 );
  if ( tol(wt,65.38,0.1)==1 ) return( 30 );
  if ( tol(wt,69.72,0.1)==1 ) return( 31 );
  if ( tol(wt,72.59,0.1)==1 ) return( 32 );
  if ( tol(wt,74.92,0.1)==1 ) return( 33 );
  if ( tol(wt,78.96,0.1)==1 ) return( 34 );
  if ( tol(wt,79.91,0.1)==1 ) return( 35 );
  if ( tol(wt,83.80,0.1)==1 ) return( 36 );
  if ( tol(wt,85.47,0.1)==1 ) return( 37 );
  if ( tol(wt,87.62,0.1)==1 ) return( 38 );
  if ( tol(wt,88.91,0.1)==1 ) return( 39 );
  if ( tol(wt,91.22,0.1)==1 ) return( 40 );
  if ( tol(wt,92.91,0.1)==1 ) return( 41 );
  if ( tol(wt,95.94,0.1)==1 ) return( 42 );
  if ( tol(wt,98.91,0.1)==1 ) return( 43 );
  if ( tol(wt,101.07,0.1)==1 ) return( 44 );
  if ( tol(wt,102.90,0.1)==1 ) return( 45 );
  if ( tol(wt,106.40,0.1)==1 ) return( 46 );
  if ( tol(wt,107.87,0.1)==1 ) return( 47 );
  if ( tol(wt,112.40,0.1)==1 ) return( 48 );
  if ( tol(wt,114.82,0.1)==1 ) return( 49 );
  if ( tol(wt,118.69,0.1)==1 ) return( 50 );
  if ( tol(wt,121.75,0.1)==1 ) return( 51 );
  if ( tol(wt,127.60,0.1)==1 ) return( 52 );
  if ( tol(wt,126.90,0.1)==1 ) return( 53 );
  if ( tol(wt,131.30,0.1)==1 ) return( 54 );
  if ( tol(wt,132.91,0.1)==1 ) return( 55 );
  if ( tol(wt,137.34,0.1)==1 ) return( 56 );
  if ( tol(wt,138.91,0.1)==1 ) return( 57 );
  if ( tol(wt,140.12,0.1)==1 ) return( 58 );
  if ( tol(wt,140.91,0.1)==1 ) return( 59 );
  if ( tol(wt,144.24,0.1)==1 ) return( 60 );
  if ( tol(wt,145.00,0.1)==1 ) return( 61 );
  if ( tol(wt,150.35,0.1)==1 ) return( 62 );
  if ( tol(wt,151.96,0.1)==1 ) return( 63 );
  if ( tol(wt,157.25,0.1)==1 ) return( 64 );
  if ( tol(wt,158.92,0.1)==1 ) return( 65 );
  if ( tol(wt,162.50,0.1)==1 ) return( 66 );
  if ( tol(wt,164.93,0.1)==1 ) return( 67 );
  if ( tol(wt,167.26,0.1)==1 ) return( 68 );
  if ( tol(wt,168.93,0.1)==1 ) return( 69 );
  if ( tol(wt,173.04,0.1)==1 ) return( 70 );
  if ( tol(wt,174.97,0.1)==1 ) return( 71 );
  if ( tol(wt,178.49,0.1)==1 ) return( 72 );
  if ( tol(wt,180.95,0.1)==1 ) return( 73 );
  if ( tol(wt,183.85,0.1)==1 ) return( 74 );
  if ( tol(wt,186.20,0.1)==1 ) return( 75 );
  if ( tol(wt,190.20,0.1)==1 ) return( 76 );
  if ( tol(wt,192.22,0.1)==1 ) return( 77 );
  if ( tol(wt,195.09,0.1)==1 ) return( 78 );
  if ( tol(wt,196.97,0.1)==1 ) return( 79 );
  if ( tol(wt,200.59,0.1)==1 ) return( 80 );
  if ( tol(wt,204.37,0.1)==1 ) return( 81 );
  if ( tol(wt,207.20,0.1)==1 ) return( 82 );
  if ( tol(wt,208.98,0.1)==1 ) return( 83 );
  if ( tol(wt,210.00,0.1)==1 ) return( 84 );
  if ( tol(wt,210.00,0.1)==1 ) return( 85 );
  if ( tol(wt,222.00,0.1)==1 ) return( 86 );
  if ( tol(wt,223.00,0.1)==1 ) return( 87 );
  if ( tol(wt,226.02,0.1)==1 ) return( 88 );
  if ( tol(wt,227.00,0.1)==1 ) return( 89 );
  if ( tol(wt,232.04,0.1)==1 ) return( 90 );
  if ( tol(wt,231.00,0.1)==1 ) return( 91 );
  if ( tol(wt,238.03,0.1)==1 ) return( 92 );
  if ( tol(wt,237.00,0.1)==1 ) return( 93 );
  if ( tol(wt,244.00,0.1)==1 ) return( 94 );
  if ( tol(wt,243.00,0.1)==1 ) return( 95 );
  if ( tol(wt,"Cm",0.1)==1 ) return( 96 );
  if ( tol(wt,"Bk",0.1)==1 ) return( 97 );
  if ( tol(wt,"Cf",0.1)==1 ) return( 98 );
  if ( tol(wt,"Es",0.1)==1 ) return( 99 );
  if ( tol(wt,"Fm",0.1)==1 ) return( 100 );
  if ( tol(wt,"Md",0.1)==1 ) return( 101 );
  if ( tol(wt,"No",0.1)==1 ) return( 102 );
  if ( tol(wt,"Lr",0.1)==1 ) return( 103 );
}

function print_lats (a,b,c,alph,beta,gamm) {
  if ( debug>0 ) printf("a b c = %20.14f%20.14f%20.14f\n",a,b,c);
  if ( debug>0 ) printf("angles = %20.14f%20.14f%20.14f\n",alph,beta,gamm);
  alph *= D2R;
  beta *= D2R;
  gamm *= D2R;
  if ( debug>0 ) printf("angrad = %20.14f%20.14f%20.14f\n",alph,beta,gamm);

  ax = a;

  bx = b * cos(gamm);
  by = b * sin(gamm);

  cx = c * cos(beta);
  cy = ( 1.0 / ( b * sin(gamm) ) ) * ( c * b * cos(alph) - c * b * cos(beta) * cos(gamm) );
  cz = sqrt( c*c - cx*cx - cy*cy );

  printf("1.0 lattice constant\n");
  printf("%20.14f%20.14f%20.14f\n",a,0.0,0.0);
  printf("%20.14f%20.14f%20.14f\n",bx,by,0.0);
  printf("%20.14f%20.14f%20.14f\n",cx,cy,cz);

}

function direct_to_cart (at,i,rb1,rb2,rb3,cart) {

    cart[3];
    
    cart[0] = at[3*i+0]*rb1[1] + at[3*i+1]*rb2[1] + at[3*i+2]*rb3[1];
    cart[1] = at[3*i+0]*rb1[2] + at[3*i+1]*rb2[2] + at[3*i+2]*rb3[2];
    cart[2] = at[3*i+0]*rb1[3] + at[3*i+1]*rb2[3] + at[3*i+2]*rb3[3];

}


# tag_value
# REQUIRED: format must be tag=value with NO SPACES
# option= 'tag_and_value' returns "tag=value"
# option= 'value_only' returns "value"
function tag_value(line,tag,option){
    # must be of the form tag=val with no separation between them
    nfields=split(line,tmp);
    for (i=1;i<nfields+1;i++){
        istag=split(tmp[i],tmp2,"=");
        if ( istag>0 ){
            for (j=1;j<istag+1;j++){
                if ( tolower(tmp2[j]) ~ tag ) {
		    if ( option=="tag_and_value" ){
			printf("%s=%s  ", tmp2[j], tmp2[j+1]);
		    } else if ( option=="value_only" ) {
			printf("%s  ", tmp2[j+1]);
		    } else {
			printf("tag_value option must be 'tag_and_value' or 'value_only'\n");
		    }
		}
            }
        }
    }
}
