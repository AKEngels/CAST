#include "coords.h"
#include <vector>
#include <string>
#include <utility> 
#include "configuration.h"
#include "global.h"
#include "energy.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "scon_matrix.h"
#include "tinker_parameters.h"
#include "tinker_refine.h"

class reaccoords
{
	public:
		
		
	coords::Coordinates *cPtr;
	coords::Coordinates const * coords;
	
	reaccoords(void);
	~reaccoords(void);
	reaccoords(coords::Coordinates *c,coords::Coordinates const* );
			
	
	static ::tinker::parameter::parameters tp;
    ::tinker::parameter::parameters cparams;
    ::tinker::refine::refined refined;


	 void bonds (void);
	 void angles (void);

	 void bonds_alteration(void);
	 void angles_alteration(void);
	 coords::Representation_3D  rmsd_align(coords::Coordinates & in);

	 std::vector <double> bonds_value,angles_value;
	 std::vector < std::vector <double> > bonds_diff,angles_diff;
	 coords::Representation_3D last_structure;
	private:





	




};