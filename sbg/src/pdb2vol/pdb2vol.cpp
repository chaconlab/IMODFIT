/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Thu Jul 22 16:15:56 CEST 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef MPI_INCLUDED
#include <mpi.h>
#endif

#include <stdlib.h>
#include <math.h>
//#include <iostream.h>


#include <libvolume/include/floatops.h>
#include <libpdb/include/Macromolecule.h>
#include <cmdl/CmdLine.h>


using namespace TCLAP;
float *bfs;


//Arguments
	char input_pdb[200],output_vol[200];
	int filter_method,name_convention,atom_convention;
	float resolution,grid_size,sigma_factor;
	bool hydrogens;
	bool BFW; // new opt form factor weights

int main(int argc, char *argv[])
 {
	//Read parameters
	void parseOptions(int argc, char** argv);
	parseOptions(argc,argv);
	init_aminoacids(Rosseta,pdb);

    vlVolume *vol,*vol2;
	Macromolecule *in_pdb;
  	ConventionNames conv;
	Convention conv_atom;
	int pad[3];
	float r_half,sigm;

	//Check if relation between resolution and grid size is correct
	sigm=resolution/2.0;  // sigma
    r_half=sigm*sqrt(log(2.0))/sqrt(3.0/2.0); // r_half
   	if ( r_half/grid_size < 1.0 )  {
    		fprintf(stdout," Increase the sampling (>=%f) and try again (increase resolution or reduce voxel size)\n>",r_half);
    		exit(0);
  	}

   	//Read macromolecule
 	in_pdb= new Macromolecule("pdb");
  	in_pdb->readPDB(input_pdb);

  	fprintf(stdout,"\npdb2map> Read pdb: %s\n",input_pdb);
  	if (!hydrogens) {
  		fprintf(stdout,"pdb2map> Removed H\n");
  		in_pdb->delete_hydrogens();
  	}

	in_pdb->geoBox();
	in_pdb->getPtrAtomPropertyBFS(&bfs);

	printf("pdb2map> Number of read proteins: %d\n",in_pdb->getLimit());
    printf("pdb2map> Number of read atoms: %d\n",in_pdb->get_num_atoms());

   //Choose filter method
   switch (filter_method)
   {
    case 0:
    	// Gaussian Expansion
    	printf("pdb2map> Gaussian expansion\n");
    	vol2=in_pdb->pdb2map_real(resolution,grid_size);
    	vol2 = FOPS::crop( vol2, 1E-12 );
    	printf("pdb2map> Electronic density: %f\n",FOPS::calc_total(vol2));
    	printf("pdb2map> Volume size: %d %d %d\n",vol2->dim().x(),vol2->dim().y(),vol2->dim().z());

    break;
    case 1:
      // Fourier filter
      printf("pdb2map> Fourier filter\n");
      int pad[3];
      //padding
      if (BFW)
      vol=in_pdb->fillVolumeBFS( grid_size, bfs );
      else
      vol=in_pdb->fillVolumeNE(grid_size);

      pad[0]= (int)( (vol->dim().x()+resolution/grid_size) /4.0);
      pad[1]= (int)( (vol->dim().y()+resolution/grid_size) /4.0);
      pad[2]= (int)( (vol->dim().z()+resolution/grid_size) /4.0);
      vol2=FOPS::padVolume(vol,vlDim(pad[0],pad[1],pad[2]),false);

	  FOPS::GaussFilter(vol2, resolution);

	  printf("pdb2map> Electronic density: %f\n",FOPS::calc_total(vol2));
      printf("pdb2map> Volume size: %d %d %d\n",vol2->dim().x(),vol2->dim().y(),vol2->dim().z());

    break;
    case 2:
     // Kernel filter
     cout<<"pdb2map> Kernel filter"<<endl;
     if (BFW)
     vol=in_pdb->fillVolumeBFS( grid_size, bfs );
     else
     vol=in_pdb->fillVolumeNE(grid_size);

     float *kernel;
     int dim_vox;
     cout<<"pdb2map> Electronic Density "<<FOPS::calc_total(vol)<<endl;
     //kernel creation and padding
     FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,vol->units().x(),resolution, sigma_factor);

     vol2 = FOPS::convoluteK(vol,kernel,dim_vox);
     FOPS::crop( vol2, 1E-12);
     // Saving files
     printf("pdb2map> Electronic Density: %f\n",FOPS::calc_total(vol2));
     printf("pdb2map> Volume size: %d %d %d\n", vol2->dim().x(),vol2->dim().y(),vol2->dim().z());
     vlPoint3f aux;
     vol2->getPosition(&aux);
     printf("pdb2map> Volume origin %f %f %f\n", aux.x(), aux.y(), aux.z());
     float cx_hi, cy_hi, cz_hi;
     FOPS::center_masses( vol2, & cx_hi, & cy_hi, & cz_hi );
     printf("pdb2map> Mass center %f %f %f\n",cx_hi,cy_hi,cz_hi);
    break;
  }

   //write file
   FOPS::writeFile(vol2,argv[2]);
   printf("pdb2map> Saving file\n");
   printf("pdb2map> Conversion realized successfully.\n");
   printf("pdb2map> \n");

   	return EXIT_SUCCESS;

 }


void parseOptions(int argc, char** argv)
{
        string temp;
       CmdLine cmd(argv[0],"   ", "1.01" );

        try {



        //
        // Define required arguments no labeled
        //


        UnlabeledValueArg<string> pdb("pdb","pdb input file","default","pdb");
        cmd.add( pdb );


        UnlabeledValueArg<string> ovol("map","volume output file ","default","map");
        cmd.add( ovol );


        UnlabeledValueArg<float> res("res","Resolution",15,"resolution");
        cmd.add( res );

        UnlabeledValueArg<float> grids("grid_size", " grid or voxel size volume",2.5,"grid size");
        cmd.add( grids );

        //
        // Define labeled arguments
        //


        ValueArg<int> nc("","nc", "aa name convention of the aminoacid atoms ( iupac(0), pdb(1))",false,0,"int");
        cmd.add( nc );

        ValueArg<int> ac("","ac", "atom convention name types ( Rosseta(0), ICM(1))",false,0,"int");
        cmd.add( ac );

        ValueArg<float> sf("","sf", "sigma factor(by default: 3)",false,3,"float");
        cmd.add( sf);

        ValueArg<bool> hyd("","hyd", "use the hydrogen atoms (by default: false)",false,false,"bool");
        cmd.add( hyd );

	    ValueArg<int> fm("f","fm", "filter resolution method:\n 0-Gaussian expansion 1-Fourier filter 2-Kernel base (default)", false,2,"int");
        cmd.add( fm );


		SwitchArg bfw("","noBFW", "disable form factor weights", true);
		cmd.add( bfw );

        // Parse the command line.
        cmd.parse(argc,argv);

        if (bfw.isSet())
        	{ BFW = false; }
        else BFW = true;


        strcpy(input_pdb,((temp=pdb.getValue()).c_str()));
        // esto no debia estar aqui...
        FILE *f;
        if ( (f=fopen(input_pdb, "r"))==NULL) {
        fprintf(stderr, "\n  Error->Cannot open file '%s'\n\n", input_pdb);
        exit(1);
        } fclose(f);


		strcpy(output_vol,((temp=ovol.getValue()).c_str()));
        resolution=res.getValue();
        grid_size=grids.getValue();
        filter_method=fm.getValue();
        sigma_factor=sf.getValue();
        name_convention=nc.getValue();
        atom_convention=ac.getValue();
        hydrogens=hyd.getValue();


        if(filter_method<0 || filter_method>1 || filter_method>2)
        	{
        	filter_method=2;
        	}

        if(name_convention<0 || name_convention>1)
        	{
        	fprintf(stderr, "\n  Error->Option incorrect '-nc %d'\n\n", name_convention);
        	exit(1);
        	}
        if(atom_convention<0 || atom_convention>1)
        	{
        	fprintf(stderr, "\n  Error->Option incorrect '-ac %d'\n\n", atom_convention);
        	exit(1);
        	}

		if(resolution<=0)
			{
        	fprintf(stderr, "\n  Error-> Increase resolution %f\n\n", resolution);
        	exit(1);
        	}
		if(grid_size<=0)
			{
        	fprintf(stderr, "\n  Error-> Increase grid_size %f\n\n", grid_size);
        	exit(1);
        	}
		if(sigma_factor<=0)
			{
        	fprintf(stderr, "\n  Error->Option incorrect '-sf %f'\n\n", sigma_factor);
        	exit(1);
        	}


        } catch ( ArgException& e )
        { fprintf(stderr,"Exception\n");
        }
        //cmd.~CmdLine();

}


