#include "reaccoord.h"


::tinker::parameter::parameters reaccoords::tp;


reaccoords::reaccoords(coords::Coordinates *cptr,coords::Coordinates const *Coord_obj)
{
	cPtr=cptr;
	coords=Coord_obj;
	
  if (!tp.valid()) tp.from_file(Config::get().get().general.paramFilename);
  std::vector<size_t> types;
  for (auto atom : (*cPtr).atoms()) scon::sorted::insert_unique(types, atom.energy_type());
  cparams = tp.contract(types);
  refined.refine(*cPtr, cparams);
}

reaccoords::~reaccoords(void)
{
}

static inline coords::Cartesian_Point center_of_geom(coords::Representation_3D const &xyz)
{
  std::size_t const N = xyz.size();
  coords::Cartesian_Point p(0);
  for (std::size_t i(0u); i < N; ++i)
  {
    p += xyz[i];
  }
  p /= coords::float_type(N);
  return p;
}


coords::Representation_3D reaccoords::rmsd_align(coords::Coordinates & in)
{
  std::size_t const N = cPtr->size();
  coords::Representation_3D  t;
  if (in.xyz().size() == N)
  {
    t.resize(N);
    coords::Cartesian_Point const cog1 = center_of_geom(cPtr->xyz());
    coords::Cartesian_Point const cog2 = center_of_geom(in.xyz());
    for (size_t l = 0; l < cPtr->size(); l++)
    {
      t.push_back((cPtr->xyz(l) - cog1) + (cog2 - cog1));
    }
  }
  return t;
}

void reaccoords::bonds(void)
{
	last_structure=cPtr->xyz();
	
 
	for (auto bond : refined.bonds())
        {
          coords::Cartesian_Point bv(cPtr->xyz(bond.atoms[0]) - cPtr->xyz(bond.atoms[1]));
		  
          double d = scon::len(bv);
		
		  bonds_value.push_back(d);
	
		}
	
}




void reaccoords::angles(void)
{

	 for (auto angle : refined.angles())
        {
          coords::Cartesian_Point 
            av1(cPtr->xyz(angle.atoms[0]) - cPtr->xyz(angle.atoms[1])),
            av2(cPtr->xyz(angle.atoms[2]) - cPtr->xyz(angle.atoms[1]));
          double const d(scon::angle_refined(av1, av2).degrees());
		      angles_value.push_back(d);
          //double const r(d*RATIOPI180);
	 }
}

void reaccoords::bonds_alteration(void)
{
	
	size_t i(0),k(0);
	bonds_diff.resize(refined.bonds().size());
	std::fstream bond1("BOND_LIST.dat",std::ios::app),bond2("BOND_MAX.dat",std::ios::app),bond3("COMPLETE_BOND_LIST.dat",std::ios::app);
	

	last_structure=cPtr->xyz();

	if(!bonds_value.empty()){
	for (auto bond : refined.bonds())
        {
			bond1<<"  index  :"<<k<<"   Atom 1   : "<<bond.atoms[0]+1<<"  Atom 2   : "<<bond.atoms[1]+1<<lineend;
			k++;
		}


	do
	{
		
		for (size_t j=0;j<refined.bonds().size();j++){

		/*std::cout<<"i  "<<i+j<<" i+ "<<i+j+refined.bonds().size()<<"    "<<bonds_value[i+j]-bonds_value[i+j+refined.bonds().size()]<<lineend;
*/
			bonds_diff[i/refined.bonds().size()].push_back(abs(bonds_value[i+j]-bonds_value[i+j+refined.bonds().size()]));
			
		}
		i+=refined.bonds().size();
	}while(i < bonds_value.size()-refined.bonds().size());


	for( i=0;i<cPtr->mult_struc_counter-1;i++)
		{

			size_t index=0;
			double max=abs(bonds_diff[i][0]);
			
			if (bonds_diff[i].empty()) continue;
			/*bond3<<i<<lineend;*/
			for (size_t j=0; j< bonds_diff[i].size();j++)
				{
					
				 
						
					if (abs(bonds_diff[i][j]) > max)
						{
							index=j;
							max = abs(bonds_diff[i][j]);
					
						}

				}
			bond2<<"INDEX  "<<index<<"max "<<max<<lineend; 
			
		}
	
	for (size_t j=0; j< bonds_diff[0].size();j++)
	{
		bond3<<j<<lineend;
	for( i=0;i<cPtr->mult_struc_counter-1;i++)
		{
	bond3<<"  "<<i<<"  "<<bonds_diff[i][j]<<lineend;
		}
	}
	}
}

void reaccoords::angles_alteration(void)
{
	std::cout<<angles_value.size()<<lineend;
	angles_diff.resize(refined.angles().size());
	size_t i(0),k(0);
	std::fstream angle1("ANGLE_LIST.dat",std::ios::app),angle2("ANGLE_MAX.dat",std::ios::app),angle3("COMPLETE_ANGLE_LIST.dat",std::ios::app);
	
	if(!angles_value.empty()){


	for (auto angle : refined.angles())
        {
			angle1<<"  index  :"<<k<<"   Atom 1   : "<<angle.atoms[0]+1<<"  Atom 2   : "<<angle.atoms[1]+1<<"  Atom 3   : "<<angle.atoms[2]+1<<lineend;
			k++;
		}

	do
	{
		
		for (size_t j=0;j<refined.angles().size();j++){

		//std::cout<<"i  "<<i+j<<" i+ "<<i+j+refined.angles().size()<<"    "<<angles_value[i+j]-angles_value[i+j+refined.angles().size()]<<lineend;
		angles_diff[i/refined.angles().size()].push_back(abs(angles_value[i+j]-angles_value[i+j+refined.angles().size()]));
		}
		i+=refined.angles().size();
	}while(i < angles_value.size()-refined.angles().size());


	for( i=0;i<cPtr->mult_struc_counter-1;i++)
		{

			size_t index=0;
			double max=abs(angles_diff[i][0]);
			
			if (angles_diff[i].empty()) continue;
			for (size_t j=0; j< angles_diff[i].size();j++)
				{
					
				 
			
					if (abs(angles_diff[i][j]) > max)
						{
							index=j;
							max = abs(angles_diff[i][j]);
					
						}

				}
			angle2<<"INDEX  "<<index<<"max "<<max<<lineend; 
		}

	for (size_t j=0; j< angles_diff[0].size();j++)
	{
		angle3<<j<<lineend;
	for( i=0;i<cPtr->mult_struc_counter-1;i++)
		{
	angle3<<"  "<<i<<"  "<<angles_diff[i][j]<<lineend;
		}
	}
	}
}