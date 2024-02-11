#!/bin/bash

#################################################
# Editable things
#################################################
# 

script_name="find_inter.sh";
version=1.4.1
ver_date="March 2017"

# version 1.4.1 Fri Mar 24 2017
#         add checksum to look for specific chemistry inters
# version 1.4.0 Thu Feb 12 10:07:34 CST 2015
#         add conversion of CHG to xyz rho and find rho at inters
# version 1.3.0 Wed Feb 11 08:51:30 CST 2015
#         when using -Z switch, print out the interstitials
#         that correspond to 'occupied' sites
# version 1.2.0 small fixes and enhancements to plotting
# version 1.1, Feb 1, 2015, add switch to ignore
#         atomic number Z in input POSCAR
#
# version 1.0, Jan 2015


#################################################
# No need to edit below this line
#################################################

FI_HOME="$HOME/src/find_interstitials"
FI="$HOME/bin/find_interstitials"

P2V4="$HOME/awkfiles/poscar_2v4.awk"
P2P="$HOME/awkfiles/poscar2powder.awk"
Pow2POS="$HOME/awkfiles/powder2poscar.awk"
CHEM="$FI_HOME/check_vertex_chem.awk"
CHG2XYZ="$HOME/awkfiles/CHG_to_xyz_rho.awk"

######################################
# defaults
calc_site_distr=0
checksum="0";
chgfile="";
dbg=0
Dstep=0.5
file=POSCAR
gen_poscar=0
halfcell=0
help=0
ignoreZ=0
ivsigmin=0.0
ivsigmax=1.0
isV5=0
make_plot=0
n_coord=4
norm=0
quiet=0
scmm=1
show_coord=0
output="inters.out"

######################################
declare SWITCH
while getopts "cC:d:D:f:GhHi:I:K:n:No:p:s:SZ:" SWITCH; do
    case $SWITCH in
    c) calc_site_distr=1; quiet=1 ;;
    C) chgfile=$OPTARG ;;
    d) dbg=$OPTARG ;;
    D) Dstep=$OPTARG ;;
    f) file=$OPTARG ;;
    G) gen_poscar=1 ;;
    h) help=1 ;;
    H) halfcell=1 ;;
    i) ivsigmin=$OPTARG ;;
    I) ivsigmax=$OPTARG ;;
    K) checksum=$OPTARG; show_coord=1 ;;
    n) n_coord=$OPTARG; output=${output%%.out}"_n"$n_coord ;;
    N) norm=1 ;;
    o) output=$OPTARG ;;
    p) make_plot=$OPTARG ;;
    s) scmm=$OPTARG ;;
    S) show_coord=1 ;;
    Z) ignoreZ=$OPTARG ;;
    esac
done

if [ $# -eq 0 ] || [ $help -eq 1 ]; then
    echo
    echo "#######################"
    echo "#      "$script_name
    echo "#######################"
    echo
    echo "Version "$version
    echo $ver_date
    echo
    echo "use: "$script_name" -f infile  [cC:d:D:f:GhHi:I:K:n:o:p:s:SZ:]"
    echo
    echo "    -f -*- infile (default POSCAR)"
    echo "    -d --- debug (debug level, default $dbg)"
    echo "    -h --- help"
    echo
    echo "    -D -*- Dstep (default $Dstep)"
    echo "           grid step size"
    echo "    -i -*- iv_sigma_min (default $ivsigmin)"
    echo "           min sigma for interstit-vertex distance variation"
    echo "    -I -*- iv_sigma_max (default $ivsigmax)"
    echo "           max sigma for interstit-vertex distance variation"
    echo "    -n -*- N_coord 4..10 (default $n_coord)"
    echo "           4=tetr 5=pent 6=oct"
    echo "    -S --- show vertex coordinates of interstitial space"
    echo "    -K -*- enter checksum to look for specific chemistry"
    echo "           of the vertex sites. Enter the sum of all vertex Z."
    echo
    echo "    -s -*- set supercell min max for search (default $scmm)"
    echo "    -H --- halfcell (restrict to 0.5 0.5 0.5)"
    echo "    -Z -*- ignore atomic number Z in input POSCAR"
    echo "           (will remove atoms Z in powder input file)"
    echo "    -C -*- give CHG file, and calculate the e-density at"
    echo "           the interstitial positions"
    echo
    echo "    -c --- calculate the distribution of sites vs iv_sig"
    echo "           output is ivsigdistr.\$output.dat"
    echo "           !!! This turns off printing to stdout of the"
    echo "           interstitials. (assuming large number here)!!!"
    echo "    -G --- generate poscar file with interstitials"
    echo "           name will be POSCAR_PLUS_INTS"
    echo "           This only cats the fract coords to the end"
    echo "           of the POSCAR file and does not alter the"
    echo "           chemistry header lines. This must be done"
    echo "           by hand."
    echo "    -o -*- output filename (default $output)"
    echo "    -p -*- make plot with gnuplot (data in inters.plt.dat)"
    echo "           arg is column number, 4=iv_ave 5=iv_sig 6=sigma"
    echo
    echo
    echo Eric Majzoub
    echo University of Missouri - St. Louis
    echo Jan 2015
    echo
    echo
    exit
fi

#############################################################
#############################################################

# What system are we running on...
system=$(uname -a | gawk '{print $1}');
if [ $system == "Darwin" ]; then
    term="x11"
else
    #    term="wxt" # broken with wxWidgets 3
    term="qt"
fi

if [ $dbg -eq 1 ]; then
    echo "--------------------------------------"
    echo "---- COMMAND LINE SWITCH SETTINGS ----"
    echo "file= "$file
    echo "checksum= "$checksum
    echo "dbg= "$dbg
    echo "Dstep= "$Dstep
    echo "gen_poscar= "$gen_poscar
    echo "halfcell= "$halfcell
    echo "help= "$help
    echo "ignoreZ= "$ignoreZ
    echo "ivsigmin= "$ivsigmin
    echo "ivsigmax= "$ivsigmax
    echo "make_plot= "$make_plot
    echo "n_coord= "$n_coord
    echo "output= "$output
    echo "scmm= "$scmm
    echo "show_coord= "$show_coord
    echo "system= "$system
    echo "term= "$term
    echo "--------------------------------------"
    echo "--------------------------------------"
fi


###################################################
# Functions
###################################################

autodetectV5 () {
    if [ ! -s $file ]; then
        isV5=0;
        echo "File "$file " is empty. Exiting."
        exit
    else
        isV5=$(cat $file | awk '(NR==6 && $0!~"[0-9]"){print 1}(NR==6 && $0~"[0-9]"){print 0}');
    fi
}


make_gplot_file () {
    
CC=$make_plot

if [ -e tmp.findinters.onlyz.cart ]; then
    SPLOT="'inters_n$n_coord.plt.dat' u 1:2:3:$CC with points pointtype 7 pointsize 2 lc palette,\
'tmp.findinters.boxpoints' u 1:2:3 with lines lt 3 lw 2,\
'tmp.findinters.onlyz.cart' u 1:2:3 with points pointtype 3 pointsize 3"
else
    SPLOT="'inters_n$n_coord.plt.dat' u 1:2:3:$CC with points pointtype 7 pointsize 2 lc palette,\
'tmp.findinters.boxpoints' u 1:2:3 with lines lt 3 lw 2"
fi

cat > gplot.p<<EOF

#set pm3d
set terminal $term
unset border
unset xtics; unset ytics; unset ztics;
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
splot $SPLOT
pause mouse keypress

EOF

}

make_control_file () {

cat > control.dat<<EOF

&condat
debug = $dbg
N_coord = $n_coord
show_coord = $show_coord
half_cell = $halfcell
scmm = $scmm
dstep = $Dstep
iv_sig_min = $ivsigmin
iv_sig_max = $ivsigmax
/

EOF

}

make_distr_awkfile () {

cat > tmp.fi.awk<<EOF

@include /home/ehm/awkfiles/awklib.awk
BEGIN{

    debug = 0;
    normalize = $norm;

    TOL = 1e-2;
    MAXDATA=10000;

    ivave[MAXDATA]=0; ivavecnt=0;
    ivsig[MAXDATA]=0; ivsigcnt=0;
    sigma[MAXDATA]=0; sigmacnt=0;

    ivsig_ave[MAXDATA]=0; # average iv_ave for each iv_sig

    nivave[MAXDATA]=0; 
    nivsig[MAXDATA]=0; 
    nsigma[MAXDATA]=0; 

    been_called=0;

}

(NF==8 && \$1=="inter"){
    
    iv_ave = \$6; iv_sig=\$7; sig=\$8;

    if (debug)
    printf("####  reading %10.4f%10.4f%10.4f\n",iv_ave,iv_sig,sig);
    
    if ( been_called == 0 ){
	if (debug) printf("setting intial values\n");
	ivave[ivavecnt++] = iv_ave; nivave[0]++;
	ivsig[ivsigcnt++] = iv_sig; nivsig[0]++;
	sigma[sigmacnt++] =    sig; nsigma[0]++;
        ivsig_ave[0] = iv_ave;
	if (debug)
	printf("ivave[0] = %10.4f  ivsig[0]= %10.4f  sigma[0]= %10.4f\n",ivave[0],ivsig[0],sigma[0]);
	been_called=1;
    } else {
	if (debug) printf("looking for repeats ...\n");
	################################
	# iv_ave
	repeat=0;
	for(i=0;i<ivavecnt;i++){
	    if (debug)
	    printf("ivave loop ivave[%3d]= %10.4f cmp %10.4f\n",i,ivave[i],iv_ave);
	    if ( tol( iv_ave , ivave[i] , TOL ) == 1 ){
		nivave[i]++;
		if (debug) printf("repeat ivave, nivave[%d]= %d\n",i,nivave[i]);
		repeat=1;
		break;
	    }
	}
	# not in the list if at this point
	if (repeat==0){
	    if (debug) printf("new ivave ... %10.4f\n",iv_ave);
            nivave[ivavecnt]++;
	    ivave[ivavecnt++] = iv_ave;
	}
        ###################################

	################################
	# iv_sig
	repeat=0;
	for(i=0;i<ivsigcnt;i++){
	    if (debug)
	    printf("ivsig loop ivsig[%3d]= %10.4f cmp %10.4f\n",i,ivsig[i],iv_sig);
	    if ( tol( iv_sig , ivsig[i] , TOL ) == 1 ){
		nivsig[i]++;
                ivsig_ave[i] += iv_ave;
		if (debug) printf("repeat ivsig, nivsig[%d]= %d\n",i,nivsig[i]);
		repeat=1;
		break;
	    }
	}
	# not in the list if at this point
	if (repeat==0){
	    if (debug) printf("new ivsig ... %10.4f\n",iv_sig);
            nivsig[ivsigcnt]++;
	    ivsig[ivsigcnt++] = iv_sig;
            ivsig_ave[ivsigcnt-1] = iv_ave;
	}
        ###################################

	################################
	# sigma
	repeat=0;
	for(i=0;i<sigmacnt;i++){
	    if (debug)
	    printf("sigma loop sigma[%3d]= %10.4f cmp %10.4f\n",i,sigma[i],sig);
	    if ( tol( sig , sigma[i] , TOL ) == 1 ){
		nsigma[i]++;
		if (debug) printf("repeat sigma, nsigma[%d]= %d\n",i,nsigma[i]);
		repeat=1;
		break;
	    }
	}
	# not in the list if at this point
	if (repeat==0){
	    if (debug) printf("new sigma ... %10.4f\n",sig);
            nsigma[sigmacnt]++;
	    sigma[sigmacnt++] = sig;
	}
        ###################################

    }

    
}

END{

    if ( debug ) {
	printf("sigma count= %d\n",sigmacnt);
	printf("ivsig count= %d\n",ivsigcnt);
	printf("ivave count= %d\n",ivavecnt);
    }


    for (i=0;i<ivavecnt;i++){if(nivave[i]>cntmax){cntmax=nivave[i]}}
    for (i=0;i<ivavecnt;i++){
        if ( normalize ) {
	  printf("%15.6f%15.6f\n", ivave[i], nivave[i]/cntmax) >> "tmp.fi.ivave.dat"
        } else {
	  printf("%15.6f%10d\n", ivave[i], nivave[i]) >> "tmp.fi.ivave.dat"
        }
    }

    for (i=0;i<ivsigcnt;i++){if(nivsig[i]>cntmax){cntmax=nivsig[i]}}
    for (i=0;i<ivsigcnt;i++){
        if ( normalize ) {
	  printf("%15.6f%15.6f\n", ivsig[i], nivsig[i]/cntmax) >> "tmp.fi.ivsig.dat"
        } else {
          printf("%15.6f%10d\n", ivsig[i], nivsig[i]) >> "tmp.fi.ivsig.dat"
        }
    }

    for (i=0;i<sigmacnt;i++){if(nsigma[i]>cntmax){cntmax=nsigma[i]}}
    for (i=0;i<sigmacnt;i++){
        if ( normalize ) {
	   printf("%15.6f%15.6f\n", sigma[i], nsigma[i]/cntmax) >> "tmp.fi.sigma.dat"
        } else {
	   printf("%15.6f%10d\n", sigma[i], nsigma[i]) >> "tmp.fi.sigma.dat"
        }
    }

    for (i=0;i<ivsigcnt;i++){
          printf("%15.6f%15.6f\n", ivsig[i], ivsig_ave[i]/nivsig[i]) >> "tmp.fi.ivsig_ave.dat"    
    }

}


EOF

}

function make_corners_awk () {
    cat > tmp.print_corners.awk<<EOF
BEGIN{
}

(NR==3){a[1]=\$1;a[2]=\$2; a[3]=\$3;}
(NR==4){b[1]=\$1;b[2]=\$2; b[3]=\$3;}
(NR==5){c[1]=\$1;c[2]=\$2; c[3]=\$3;}

END{
printf("%12.5f%12.5f%12.5f\n",0,0,0);
printf("%12.5f%12.5f%12.5f\n",a[1],a[2],a[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+b[1],a[2]+b[2],a[3]+b[3]);
printf("%12.5f%12.5f%12.5f\n",b[1],b[2],b[3]);
printf("%12.5f%12.5f%12.5f\n",0,0,0);
printf("%12.5f%12.5f%12.5f\n",c[1],c[2],c[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+c[1],a[2]+c[2],a[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+b[1]+c[1],a[2]+b[2]+c[2],a[3]+b[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",b[1]+c[1],b[2]+c[2],b[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",c[1],c[2],c[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+c[1],a[2]+c[2],a[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",a[1],a[2],a[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+b[1],a[2]+b[2],a[3]+b[3]);
printf("%12.5f%12.5f%12.5f\n",a[1]+b[1]+c[1],a[2]+b[2]+c[2],a[3]+b[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",b[1]+c[1],b[2]+c[2],b[3]+c[3]);
printf("%12.5f%12.5f%12.5f\n",b[1],b[2],b[3]);
}

EOF
}

function make_inter_correlate (){
cat > inter_corr.awk <<EOF

BEGIN{
    # This script will take the output of find_interstitials.f90
    # and command line input fx,fy,fz (a particular fractional
    # coordinate) and find the nearest interstitial position
    # to fx, fy, fz. It uses the basis to calculate the distance
    # in angstroms (dij in output).

    # set on command line
    #fx; fy; fz;
    #debug=0;
    ic=1;
    dij_min=1e9;
    ib[10]=0; icnt=1;
    dij[10000]=0;
}

(NF==4 && \$1=="IB:"){
    ib[icnt]=\$2; ib[icnt+1]=\$3; ib[icnt+2]=\$4;
    icnt=icnt+3;
}

(NF==5 && \$1=="frac"){
    inter[ic]=\$2; x[ic]=\$3; y[ic]=\$4; z[ic]=\$5;
    X[ic] = x[ic] * ( ib[1] + ib[4] + ib[7] );
    Y[ic] = y[ic] * ( ib[2] + ib[5] + ib[8] );
    Z[ic] = z[ic] * ( ib[3] + ib[6] + ib[9] );
    ic++;
}

END{
    if ( debug>0 ) {
      printf("input basis:\n");
      printf("%10.5f%10.5f%10.5f\n",ib[1],ib[2],ib[3]);
      printf("%10.5f%10.5f%10.5f\n",ib[4],ib[5],ib[6]);
      printf("%10.5f%10.5f%10.5f\n",ib[7],ib[8],ib[9]);
    }

    FX = fx * ( ib[1] + ib[4] + ib[7] );
    FY = fy * ( ib[2] + ib[5] + ib[8] );
    FZ = fz * ( ib[3] + ib[6] + ib[9] );

    for(i=1;i<ic;i++){
	dij[i] = sqrt((FX-X[i])^2+(FY-Y[i])^2+(FZ-Z[i])^2);
        if ( debug>0 ) printf("i= %4d   dij= %8.4f  %8.4f%8.4f%8.4f  %8.4f%8.4f%8.4f\n",i,dij[i],x[i],y[i],z[i],fx,fy,fz);
	if ( dij[i] < dij_min ) {
	    dij_min=dij[i]; internum=inter[i];
	}
    }
    printf("%d  %f",internum,dij[internum]);
    #printf("closest to interstitial %d  dij= %f\n",internum,dij[internum]);
}

EOF
}

function make_chg_corr_awk () {
cat > chg_corr.awk<<EOF
BEGIN{

dijmin2=1e9;
rhomin=1e9;
never_found=1;

dtol=0.25; # frac coords
cxmd=cx-dtol;
cxpd=cx+dtol;
cymd=cy-dtol;
cypd=cy+dtol;
czmd=cz-dtol;
czpd=cz+dtol;

}
(NF==4){
gx=\$1; gy=\$2; gz=\$3; rho=\$4;

if ( gx<cxmd ) next;
if ( gx>cxpd ) next;
if ( gy<cymd ) next;
if ( gy>cypd ) next;
if ( gz<czmd ) next;
if ( gz>czpd ) next;

if ( debug ) printf("gx gy gz cx cy cz %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",gx,gy,gz,cx,cy,cz);
dij2 = (cx-gx)^2 + (cy-gy)^2 + (cz-gz)^2;
  if ( dij2 < dijmin2 ){
    never_found=0;
    dijmin2=dij2; rhomin=rho;
  }
}
END{
 # There is a bug here. If the cart coords would give something
 # outside of 0..1, then the grid points (all in 0..1) may not
 # find this point. Fixable, but rare, so letting it go for now.
  if ( never_found ) printf("%15.8f%15.8f\n",-1,-1);
  else printf("%15.8g%15.8g\n",sqrt(dijmin2),rhomin);
}
EOF
}

clean_up () {

    rm -f fi_in.dat tmp.fi.*;
    rm -f control.dat;
    rm -f tempfile.poscar_cnvrt.POSv4;
    rm -f tmp.findinters.onlyz.pos;
    rm -f tmp.findinters.onlyz.cart;
    rm -f tmp.findinters.onlyz.frac;
    rm -f tmp.findinters.boxpoints;
    rm -f tmp.inter_corr.dat;
    rm -f tmp.print_corners.awk;
    rm -f gplot.p;
    rm -f inters_n*.plt.dat;
    rm -f tmp.findinters.chgxyz.dat;
    rm -f tmp.chgdat.*;
    rm -f chg_corr.awk;
    rm -f inter_corr.awk;
}

##################################
# BEGIN SCRIPT
##################################
autodetectV5;
OPTS_P2P="-v ignoreZ=$ignoreZ"
if [ $isV5 -eq 1 ]; then
    cat $file | igawk -f $P2V4  | igawk -f $P2P $OPTS_P2P > fi_in.dat;
else
    cat $file | igawk -f $P2P $OPTS_P2P > fi_in.dat;
fi

make_control_file;

echo "Running find_interstitials..."
if [ $quiet -eq 1 ]; then
    $FI > $output;
else
    $FI | tee $output;
fi

echo
echo " +++ fortran run finished +++"
echo
echo " === processed output below ==="
echo

cat $output | gawk '($1=="min"){print $0}'
cat $output | gawk '($1=="sigma"){print $0}'
cat $output | gawk '($1=="d_ave"){print $0}'
cat $output | gawk '($1=="iv_ave"){print $0}'
cat $output | gawk 'BEGIN{iv_sig_min=1e9;iv_sig_max=0;}($1=="inter"){if($7>iv_sig_max){iv_sig_max=$7}if($7<iv_sig_min){iv_sig_min=$7}}END{printf("     iv_sig     %10.6f%15.6f\n",iv_sig_min,iv_sig_max);}'

echo
echo "  use   -i iv_sig_min    and    -I iv_sig_max"
echo

if [ $calc_site_distr -eq 1 ]; then
    echo "Calculating site distribution..."
    make_distr_awkfile;
    cat $output | igawk -f ./tmp.fi.awk;
    mv -f tmp.fi.ivave.dat $output.ivave.dat;
    mv -f tmp.fi.ivsig.dat $output.ivsig.dat;
    mv -f tmp.fi.sigma.dat $output.sigma.dat;
    mv -f tmp.fi.ivsig_ave.dat $output.ivsig_vs_ave_iv_ave.dat;
    echo 
fi

if [ $gen_poscar -eq 1 ]; then
    cat $file > POSCAR_PLUS_INTS;
    cat $output | grep frac | gawk '{print $3,$4,$5}' >> POSCAR_PLUS_INTS;
fi

#################################################################
if [ $ignoreZ -gt 0 ]; then
    if [ $isV5 -eq 1 ]; then
	cat $file | igawk -f $P2V4 | igawk -f $P2P -v onlyZ=$ignoreZ | igawk -f $Pow2POS -v PFT=3 > tmp.findinters.onlyz.pos;
    else
	cat $file | igawk -f $P2P -v onlyZ=$ignoreZ | igawk -f $Pow2POS -v PFT=3 > tmp.findinters.onlyz.pos;
    fi
    
    cat tmp.findinters.onlyz.pos | gawk --source '(NR>6 && NF==3){print $0}' > tmp.findinters.onlyz.frac;
    poscar_cnvrt -f tmp.findinters.onlyz.pos -C | gawk --source '(NR>6 && NF==3){print $0}' > tmp.findinters.onlyz.cart;

    # find the occupied interstitials
    echo "       -------------- Occupied Interstitials Sites Z= "$ignoreZ" -------------------"
    echo
    echo "                        x        y        z  iv_ave  iv_sig   sigma      d-i"
    make_inter_correlate;
    cat tmp.findinters.onlyz.frac | while read ix iy iz; do
	if [ $dbg -gt 0 ]; then
	    echo "fx= "$ix "  fy= "$iy "  fz= "$iz
	fi
	cat $output | igawk -f ./inter_corr.awk -v fx=$ix -v fy=$iy -v fz=$iz > tmp.inter_corr.dat;
	inum=$(cat tmp.inter_corr.dat | gawk '{print $1}')
	dij=$(cat tmp.inter_corr.dat | gawk '{print $2}')
	cat $output | gawk --source '($1=="inter" && $2==INUM){printf("%s%9.4f\n",$0, DIJ);}' -v INUM=$inum -v DIJ=$dij
    done
    echo
    echo "    note: The d-i (distance to interstitial site) values should be near zero."
    echo "          If not, the coordination number is likely incorrect. "
    echo
fi


#################################################################
# Get the electron density at every inter if CHG file is present
if [ "$chgfile" != "" ]; then
    echo "Converting CHG to xyz_rho file..."
    echo
    make_chg_corr_awk;
    cat $chgfile | igawk -f $CHG2XYZ > tmp.findinters.chgxyz.dat;
    TMP="       ------------------ Charge Density at Interstitial Positions ---------------"
    echo "$TMP"; echo "$TMP"  > rho_n$n_coord.dat
    echo
    TMP="              nn        x        y        z  iv-ave  iv-sig   sigma   grid-int      rho"
    echo "$TMP"; echo "$TMP" >> rho_n$n_coord.dat

    cat $output | gawk --source '($1=="inter"){print $2,$3,$4,$5}' | while read nn xx yy zz; do
	cat tmp.findinters.chgxyz.dat | igawk -f ./chg_corr.awk -v cx=$xx -v cy=$yy -v cz=$zz > tmp.chgdat.$nn;
	DIJ=$(cat tmp.chgdat.$nn | gawk '{print $1}')
	RHO=$(cat tmp.chgdat.$nn | gawk '{print $2}')
	cat $output | gawk --source '($1=="inter" && $2==NN){printf("%s%10.4f%10.4f\n",$0,dd,rr);}' -v NN=$nn -v dd=$DIJ -v rr=$RHO
    done | tee -a rho_n$n_coord.dat
    echo
    echo "    note: If grid-int and rho are -1, then the cart coords were out of the"
    echo "          unit cell (frac coords were outside of 0..1). This is fixable."
    echo
fi

#################################################################

if [ $checksum -gt 0 ]; then
    echo "Running checksum on the interstitials..."
    cat inters_n$n_coord | gawk -f $CHEM -v checksum=$checksum\
	| grep "####" | sort -n -k 6
fi

if [ $make_plot -gt 0 ]; then
    
    make_gplot_file;
    make_corners_awk;
    cat $file | igawk -f ./tmp.print_corners.awk > tmp.findinters.boxpoints;
    
    echo "#        x         y         z    iv_ave    iv_sig     sigma" > inters_n$n_coord.plt.dat;
    cat $output | gawk '($1=="inter"){printf("%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",$3,$4,$5,$6,$7,$8);}' >> inters_n$n_coord.plt.dat
    
    if [ $system == "Darwin" ]; then
	echo "on Darwin system"
	echo "load 'gplot.p'; pause mouse keypress" | gnuplot -persist
    else
	gnuplot  -persist  gplot.p 
    fi
fi

echo " +++ Done +++"

if [ $dbg -eq 0 ]; then
    echo "cleaning up..."
    echo
    clean_up;
fi


exit
