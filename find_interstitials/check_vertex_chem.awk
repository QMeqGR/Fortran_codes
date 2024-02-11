BEGIN{

    # SET ON COMMAND LINE
    #checksum - the sum of all Z you are looking for
    read_from_lines=1;
    vnum=0;

}

($2=="N_coord"){N_coord=$3;}
($1=="inter"){
    printf("inter %3d ivave= %f  ivsig= %f ",$2,$6,$7);
}
($1=="from"){
    Z[vnum++]=$7;
    if (vnum == N_coord) {
	vnum=0; tot=0;
	for(i=0; i<N_coord; i++){
	    tot += Z[i];
	    printf("%3d ",Z[i]);
	}
	printf(" chksum= %d ",tot);
	if ( tot == checksum ){
	    printf("####\n");
	} else {
	    printf("\n");
	}
    }
}

END{

}
