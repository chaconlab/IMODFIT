/*
  * voltool.cpp
 *
 *  Created on: Mar 19, 2018
 *      Author: pablo
 */

#define prog "voltool" // Program name
#define VERSION "1.1" // Version code

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <libtools/include/timer.h>
#include <libvolume/include/vlvolume.h>
#include <libvolume/include/MarchingCubes.h>
#include <libvolume/include/floatops.h>
#include <libpdb/include/Macromolecule.h>

#include <cmdl/CmdLine.h>


using namespace TCLAP;

char input_map[256],output_map[256],extension[200],*point, in_maskpdb[256];
float euler[3],offset[3];
int euler_convention;
bool in_pad;
bool in_endian;
bool in_origin;
bool in_mask;
bool in_cut=false;
float in_surf=10;

float in_cutoff,in_vox,in_res, voxelsize, in_res_mask;

int main(int argc, char *argv[])
{
	//Read parameters
	void parseOptions(int argc, char** argv);
	parseOptions(argc,argv);

	vlPoint3f pos, pos2;

	timer t;
	float ax,ay,az,auxX,auxY,auxZ;
	float centervol[3];
	float voxelSize=-1;
	bool crop=false;
	int i;
	char aux;

	//
	// Reading
	//
	vlVolume *vl, *vol, *vol_dump;
	fprintf( stdout, "%s> Reading %s\n", prog, input_map);

	if(voxelSize>0)
		vl=FOPS::readFile(input_map,in_vox);
	else
		vl=FOPS::readFile(input_map);
	if(vl==NULL)
	{
		printf("Error: unknown inputfile type\n");
		return 1;
	}
	fprintf( stdout, "%s> Dimensions %dx%dx%d  voxsize %f\n",prog, vl->dim().x(),
			vl->dim().y(), vl->dim().z(), vl->units().x()  );

	//
	// Density threshold & Crop
	//
	if (in_cut) {
		fprintf( stdout, "%s> Density threshold & Crop %f\n", prog, in_cutoff);
		FOPS::threshold( vl, in_cutoff);
		vl=FOPS::crop( vl, 1E-15 ,true);
		fprintf( stdout, "%s> Resize %dx%dx%d\n",prog, vl->dim().x(),
				vl->dim().y(), vl->dim().z()  );
		fprintf( stdout, "%s> Density below the cutoff  %.3f \n", prog, FOPS::calc_total(vl));
	}
	vl->getPosition(&pos);






	t.restart();

	//
	// Volumen Centering
	//
	if(in_origin)
	{
		FOPS::center_vol(vl,centervol);
		vlPoint3f aux_pos;
		aux_pos.x(centervol[0]*vl->units().x());
		aux_pos.y(centervol[1]*vl->units().y());
		aux_pos.z(centervol[2]*vl->units().z());
		//	pos2.x(-aux_pos.x()+offset[0]);
		//	pos2.y(-aux_pos.y()+offset[1]);
		//	pos2.z(-aux_pos.z()+offset[2]);
		//vol->setPosition(pos2);
		vl->setPosition(aux_pos);
		fprintf( stdout, "%s> center %.3fx%.3fx%.3f\n", prog, aux_pos.x(),aux_pos.y(),aux_pos.z());
	}





	//
	// Rotate
	//
	if(euler[0]!=0 || euler[1]!=0 || euler[2]!=0) {

		fprintf( stdout, "%s> Rotate %.3fx%.3fx%.3f\n", prog, euler[0],euler[1],euler[2]);

		ax=euler[0]*3.1415927/180.0;
		ay=euler[1]*3.1415927/180.0;
		az=euler[2]*3.1415927/180.0;

		// resize square for rotate
		//float max=FOPS::max_lenght(vl);
		vol=FOPS::padVolume(vl,vlDim((vl->dim().x()/2),(vl->dim().y()/2),vl->dim().z()/2));
		vol->getPosition(&pos);
		vol_dump = new vlVolume(vol);
		vol_dump->setPosition(pos);
		vol_dump->clear(0);

		if(ax==0.0)
			ax=0.0000001;
		if(ay==0.0)
			ay=0.00000001;
		if(az==0.0)
			az=0.0000009;

		FOPS::rotate_interp(vol, vol_dump, ax,ay,az, 0);
		// FOPS::rotate(vol,ax,ay,az, 0,0,0,0);

		vl=vol_dump;
		vl=FOPS::crop( vl, 1E-15 ,true);

	}

	//
	// translate
	//
	if(offset[0]!=0 || offset[1]!=0 || offset[2]!=0)
	{
		vl->getPosition(&pos);
		fprintf( stdout, "%s> Translate %.3fx%.3fx%.3f\n", prog, offset[0],offset[1],offset[2]);
		pos2.x(pos.x()+offset[0]);
		pos2.y(pos.y()+offset[1]);
		pos2.z(pos.z()+offset[2]);
		vl->setPosition(pos2);
	}



	//Kernel filter
	if(in_res>0)
	{
		float *kernel;
		int dim_vox;
		float dens0,dens1;
		FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,vl->units().x(),in_res, 3.0);
		// vlVolume *padded = FOPS::padVolume( vl , vlDim( (dim_vox-1)/2, (dim_vox-1)/2,  (dim_vox-1)/2));
		// vlVolume *padded2 = FOPS::convoluteK(padded,kernel,dim_vox);
		// FOPS::projectVolume(padded2,vl
		fprintf( stdout, "%s> Lowering resolution to  %.3f (%d)\n", prog, in_res, dim_vox);
		vl->getPosition(&pos);
		//fprintf( stdout, "%s> pos 1 %.3fx%.3fx%.3f\n", prog, pos.x(),pos.y(), pos.z());

		vol=FOPS::padVolume(vl,vlDim(vl->dim().x()/2, vl->dim().y()/2, vl->dim().z()/2));
//		fprintf( stdout, "%s>  Inital dimensions %dx%dx%d  voxsize %f\n",prog, vl->dim().x(),
//					vl->dim().y(), vl->dim().z(), vl->units().x() );
//		fprintf( stdout, "%s>  Final dimensions %dx%dx%d  voxsize %f\n",prog, vol->dim().x(),
//					vol->dim().y(), vol->dim().z(), vol->units().x() );
//		fprintf( stdout, "%s> pos 2 %.3fx%.3fx%.3f\n", prog, pos.x(),pos.y(), pos.z());

//		FOPS::threshold( vl, in_cutoff);
//		dens0=FOPS::calc_total_dens(vl);


		vol->getPosition(&pos);
		FOPS::GaussFilter(vol, in_res);
		float diag=sqrt( pow(vl->dim().x(),2) + pow(vl->dim().y(),2) + pow(vl->dim().z(),2))/2.0;
		//FOPS::GaussBandFilter(vol, in_res, diag );
		//FOPS::ButterFilter(vol, in_res, 2, 0);





		//fprintf( stdout, "%s> pos 3 %.3fx%.3fx%.3f\n", prog, pos.x(),pos.y(), pos.z());

		vl->setPosition(pos);
		vl->clear(0);
		vl=FOPS::resize(vol,vlDim(vl->dim().x(),vl->dim().y(),vl->dim().z()));
		vl->getPosition(&pos);
		//fprintf( stdout, "%s> pos 4 %.3fx%.3fx%.3f\n", prog, pos.x(),pos.y(), pos.z());
		FOPS::threshold( vl, 1E-15 );

		vl=FOPS::crop( vl, 1E-15 ,true);

		fprintf( stdout, "%s>  Final dimensions %dx%dx%d  voxsize %f\n",prog, vl->dim().x(),
				vl->dim().y(), vl->dim().z(), vl->units().x() );

		 // fprintf( stdout, "%s> Density below the cutoff  %.3f \n", prog, FOPS::calc_total(vl));
		 // FOPS::norm_threshold( vl, in_cutoff,dens0);
		 // FOPS::threshold( vl, in_cutoff);
		 // dens0=FOPS::calc_total_dens(vl);
	     // fprintf( stdout, "%s> Density below the cutoff  %.3f \n", prog, dens0);


	}


	if(in_mask)
	{

		Macromolecule *mol;
		mol = new Macromolecule("mask");
		fprintf( stdout, "%s> Generating a mask around %s pdb\n",prog, in_maskpdb);
		mol->readPDB(in_maskpdb);

		vlVolume *map2 = new vlVolume( vl );
		map2->clear(0);
		mol->project_Real( map2 ); // projects PDB into a volume with the same size as input one

		// Computing Gaussian Kernel (filtering)
		fprintf( stdout, "%s>  Kernel filter (rmask= %f)\n",prog,in_res_mask);
		float *kernel;
		int dim_vox;
		int sigma_factor = 3;
		FOPS::compute_kernel_Gaussian(&kernel,&dim_vox,map2->units().x(),in_res_mask, sigma_factor);
		int border;
        border= (dim_vox-1)/2+in_surf;


		// pading pdb-map
		map2 = FOPS::padVolume( map2 , vlDim( border, border, border) );
		// also padding the original (initial) map (in order both maps share the same indexing)
		vl = FOPS::padVolume( vl , vlDim( border, border, border));
		vl->getPosition(&pos);

		// PDB MAP FILTRATION
		fprintf( stdout, "%s>  Map convolution (kernel size = %d)\n",prog,dim_vox);

		vlVolume *dummy = new vlVolume( map2 ); // makes copy from target map
		// dummy->clear(0);
		FOPS::convoluteK_nopad(dummy,map2,kernel,dim_vox, false);
		map2->setPosition(pos);


		S_Grid grid;
		vlVolume * mask;

		fprintf( stdout, "%s>  Surface around map R %f (surfR)\n", prog, in_surf);
		grid=FOPS::surfaceMask_grid( map2, in_surf);

		map2->clear(0);

		for(int m=0;m<map2->dim().x()+0;m++)
			for (int n=0;n<map2->dim().y()+0;n++)
				for (int k=0;k< map2->dim().z()+0;k++)
				{
					if((grid)->matrix[m][n][k].phi==0)
						map2->setVoxel(vlPoint3ui(m,n,k),(float)1.0);
				}

		float *container= (float*)vl->getVoxelVoidPtr(vlPoint3ui(0,0,0));
		float *container2= (float*)map2->getVoxelVoidPtr(vlPoint3ui(0,0,0));
		vlStep step=map2->stepping();
		int stepx=step.x();
		int stepy=step.y();
		int stepz=step.z();
		vlDim dim=map2->dim();
		int frame=(dim_vox-1)/2;
		int offset;
		float voxel,voxel2;
		for(int i=0;i<dim.x();i++)
			for(int j=0;j<dim.y();j++)
				for(int z=0;z<dim.z();z++)
				{
					offset = (i*stepx)+(j*stepy)+(z*stepz);
					voxel2= *(container2+offset); // 1st vol voxel
					if( voxel2 > 0 ) // "vol" >0 mask
						*(container2+offset) = *(container+offset); // asigning value from the original map
					else
						*(container2+offset) = 0.0; // asigning zero
				}

		map2->setPosition(pos);

		vl=map2;


		FOPS::threshold( vl, 1E-15 );

		vl=FOPS::crop( vl, 1E-15 ,true);


		// FOPS::writeFile(map2,"temp2.mrc");



	}


	//
	// Interpolating
	//

	if(voxelsize>0)
	{
		fprintf( stdout, "%s> Interpolating %.3f -> %.3f\n", prog, vl->units().x(), voxelsize);
		if(vl->units().x()< voxelsize)  {
			vl=FOPS::interpolate_map(vl, vlUnit( voxelsize, voxelsize,  voxelsize), true);
		} else {
			fprintf( stdout, "%s> Warning %.3f > %.3f doing nothing\n", prog, vl->units().x(), voxelsize);
		}
		fprintf( stdout, "%s>  Final dimensions %dx%dx%d  voxsize %f\n",prog, vl->dim().x(),
				vl->dim().y(), vl->dim().z(), voxelsize  );

	}




	//Check output extension
	point=strchr(output_map,'.');
	if(point==NULL)
	{
		std::cout<<"vol_conv> ERROR. Not recognized File extension: "<<output_map<<std::endl;
		return NULL;
	}
	i=0;
	while(point[i]!='\0')
	{
		extension[i]=point[i];
		i++;
	}
	extension[i]='\0';

	//Creates and writes the volume surface
	if(strcasecmp(extension,".vmd")==0)
	{
		MarchingCubes *mc=new MarchingCubes(vl,in_cutoff);
		mc->init_all();
		mc->run();
		mc->clean_temps();
		mc->writeVMD(argv[2]);
		mc->clean_all();
		fprintf( stdout, "%s> Write %s \n", prog, output_map);

	}
	else if(strcasecmp(extension,".wrl")==0)
	{
		MarchingCubes *mc=new MarchingCubes(vl,in_cutoff);
		mc->init_all();
		mc->run();
		mc->clean_temps();
		mc->writeVRML(argv[2]);
		mc->clean_all();
		fprintf( stdout, "%s> Write %s \n", prog, output_map);

	}

	else //Writes the volume
	{
		if(!FOPS::writeFile(vl,output_map))
		{
			fprintf( stderr, "%s> ERROR. Error: unknown outputfile type %s \n", prog, output_map);
			return 1;
		}
		fprintf( stdout, "%s> Write %s \n%s> All done \n%s>\n", prog, output_map, prog, prog);

	}

}

inline vector<string> split(const string& str, const string& delim)
		{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos-prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	}
	while (pos < str.length() && prev < str.length());
	return tokens;
		}



void parseOptions(int argc, char** argv)
{
	string temp;
	CmdLine cmd(prog,"EM volume tool", VERSION );

	try {

		//
		// Define required arguments no labeled
		//
		UnlabeledValueArg<string> map_in("map_in","map input file","default","map_in");
		cmd.add( map_in );

		UnlabeledValueArg<string> map_out("map_out","map output file","default","map_out");
		cmd.add( map_out );

		ValueArg<std::string> rotE("","rot", " Euler Rotation",false,"","Euler_1,Euler_2,Euler_3");
		cmd.add( rotE );

		ValueArg<std::string> trans("","xyz", " Translation shift (after rotation)",false,"","X,Y,Z");
		cmd.add( trans);



		// Define labeled arguments
		//

		ValueArg<float> cut("","cutoff", "Density cutoff for setting the volume surface",false,0,"float");
		cmd.add( cut );

		ValueArg<string> Pdb("p","pdbmask", "Mask pdb input file", false,"","pdb_in");
		cmd.add( Pdb );

		ValueArg<float> resM("","rmask", " Resolution of the mask pdb", false, 10,"float");
		cmd.add( resM );

		ValueArg<float> surfR("","Rmask", " Surface radius of the mask pdb", false, 10,"float");
		cmd.add( surfR );

		ValueArg<float> voxelS("s","voxel", "Voxel size of the ouput map",false,0,"float");
		cmd.add( voxelS );

		ValueArg<float> res("r","res", "Resolution to filter the map",false,-1,"float");
		cmd.add( res );

		SwitchArg pad("","pad", "activation of padding and cropping of the volume", true);
		cmd.add( pad );

		ValueArg<float> vox("","vox", "(by default: 0) voxel size (necessary if the input map is in MAT format)",false,0,"float");
		cmd.add( vox );

		SwitchArg endian( "e", "endian", "change endian (only for CCP4)", false );
		cmd.add( endian );

		SwitchArg origin( "o", "origin", "movement from origin", false );
		cmd.add( origin );

		ValueArg<int> Ec("","Ec", "(by default: 0) Convention of the Euler angles (0=x-convention 1=xyz-convention 2=y-convention)",false,0,"int");
		cmd.add( Ec );

		// Parse the command line.
		cmd.parse(argc,argv);

		strcpy(input_map,((temp=map_in.getValue()).c_str()));
		// esto no debia estar aqui...
		FILE *f;
		f=fopen(input_map, "r");
		if ( (f)==NULL) {
			fprintf(stderr, "\n  Error->Cannot open file '%s'\n\n", input_map);
			exit(1);
		} fclose(f);

		strcpy(output_map,((temp=map_out.getValue()).c_str()));

		if (Pdb.isSet()) {
			strcpy(in_maskpdb,((temp=Pdb.getValue()).c_str()));
			in_mask=true;
		}

		char dumpS[256]; // Output files basename
		vector<string>  Stemp;

		if (rotE.isSet()) {
			strcpy(dumpS,((temp=rotE.getValue()).c_str()));
			Stemp = split(dumpS, ",");
			euler[0]=strtof(Stemp[0].c_str(),0);
			euler[1]=strtof(Stemp[1].c_str(),0);
			euler[2]=strtof(Stemp[2].c_str(),0);
		} else
			for ( int i = 0; i  < 3; i++ ) euler[i] = 0 ;

		if (trans.isSet()) {
			strcpy(dumpS,((temp=trans.getValue()).c_str()));
			Stemp = split(dumpS, ",");
			offset[0]=strtof(Stemp[0].c_str(),0);
			offset[1]=strtof(Stemp[1].c_str(),0);
			offset[2]=strtof(Stemp[2].c_str(),0);
		} else
			for ( int i = 0; i  < 3; i++ ) offset[i] =0;


		euler_convention=Ec.getValue();
		if (pad.isSet())
			in_pad=true;
		else
			in_pad=false;

		in_vox=vox.getValue();
		if (cut.isSet()) in_cut=true;
		in_cutoff=cut.getValue();

		in_endian=endian.getValue();
		in_origin=origin.getValue();
		in_res=res.getValue();
		in_res_mask=resM.getValue();
		in_surf=surfR.getValue();

		voxelsize=voxelS.getValue();

		if(in_vox<0)
		{
			fprintf(stderr,"Error: incorrect vox size= %f\n", in_vox);
			exit(1);
		}

		if(in_cutoff<0)
		{
			fprintf(stderr,"Error: incorrect cutoff= %f\n", in_cutoff);
			exit(1);
		}


	} catch ( ArgException& e )
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
	}
	//cmd.~CmdLine();

}


