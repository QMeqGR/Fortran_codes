@include "/home/ehm/awkfiles/awklib.awk"
BEGIN{
    # version 2.0, 13 Feb 2015
    # fractional output option
    #
    # version 1.0, 12 Feb 2015
    # x is the fastest index , z the slowest, according to VASP website
    MAX=200;
    MAXGRID=MAX*MAX*MAX;
    direct=0;
    cart=0;
    frac_out=0;
    nat=0;
    rhocount=0; rc=0;
    rho[MAXGRID]=0;
}

( NR==2 ) { scale=$1; }
( NR==3 && NF==3 ){ a[0]=$1; a[1]=$2; a[2]=$3; }
( NR==4 && NF==3 ){ b[0]=$1; b[1]=$2; b[2]=$3; }
( NR==5 && NF==3 ){ c[0]=$1; c[1]=$2; c[2]=$3; }
( NR==6 ) {for(i=1;i<NF+1;i++){elm[i]=$i;}}
( NR==7 ) {
    for(i=1;i<NF+1;i++){
	n_elm[i]=$i;
	nat += $i;
    }
}
( NR==8 && $1~"D" ){ direct=1; }
( NR==8 && $1~"C" ){ cart=1; }
( NR==(8+nat+2) && NF=3 ){
    nx=$1; ny=$2; nz=$3;
#    printf("nx= %d  ny= %d  nz=%d\n",nx,ny,nz);
    if ( nx>MAX || ny>MAX || nz>MAX ){printf("Increase MAX in awk file. Exiting.\n"); exceedmax=1; }
}
( NR >(8+nat+2) && NF>0 ){
    if ( exceedmax==1 ) {
	printf("nx= %d  ny= %d  nz=%d   MAX= %d\n",nx,ny,nz,MAX);
	exit;
    }
    for(i=1;i<NF+1;i++){
	rho[ rhocount++ ] = $i;
    }
}

END{
    if ( exceedmax==1 ) exit;
    cv = cellvolume(a,b,c);
    # This is the matrix transpose needed
    # to find the fractional coordinates
    # ax bx cx
    # ay by cy
    # az bz cz
    M[0]=a[0];
    M[1]=b[0];
    M[2]=c[0];
    M[3]=a[1];
    M[4]=b[1];
    M[5]=c[1];
    M[6]=a[2];
    M[7]=b[2];
    M[8]=c[2];
    inverse3x3(M,Minv);

    if ( debug ) {
	printf("cell volume = %f\n",cv);
	printf("%12.6f%12.6f%12.6f\n",a[0],a[1],a[2]);
	printf("%12.6f%12.6f%12.6f\n",b[0],b[1],b[2]);
	printf("%12.6f%12.6f%12.6f\n",c[0],c[1],c[2]);
	printf("gridcount from header: %d %d %d --> %d\n",nx,ny,nz,nx*ny*nz);
	printf("read in rhocount= %d\n",rhocount);
	print_matrix_3x3(M);
	print_matrix_3x3(Minv);
    }
    for(iz=1; iz<nz+1; iz++){
	for(iy=1; iy<ny+1; iy++){
	    for(ix=1; ix<nx+1; ix++){
		x = a[0]*ix/nx + b[0]*iy/ny + c[0]*iz/nz;
		y = a[1]*ix/nx + b[1]*iy/ny + c[1]*iz/nz;
		z = a[2]*ix/nx + b[2]*iy/ny + c[2]*iz/nz;
		if ( frac_out == 1 ){
		    vin[0]=x; vin[1]=y; vin[2]=z;
		    matrix_33_vec_mult(Minv,vin,fr)
		    printf("%12.4f%12.4f%12.4f%12.4f\n",
			   fr[0],fr[1],fr[2],rho[rc++]/cv);
		}
		else {
		    printf("%12.4f%12.4f%12.4f%12.4f\n",x,y,z,rho[rc++]/cv);
		}
	    }
	}
    }
}
