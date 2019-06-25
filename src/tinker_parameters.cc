#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "tinker_parameters.h"
#include <set>
#include "filemanipulation.h"
#include "Scon/scon_utility.h"
#include <iterator>

#ifdef _MSC_VER
#pragma warning(disable: 4996)
#endif

/*

######## #### ##    ## ########  ######## ########  
##        ##  ###   ## ##     ## ##       ##     ## 
##        ##  ####  ## ##     ## ##       ##     ## 
######    ##  ## ## ## ##     ## ######   ########  
##        ##  ##  #### ##     ## ##       ##   ##   
##        ##  ##   ### ##     ## ##       ##    ##  
##       #### ##    ## ########  ######## ##     ## 

*/


tinker::parameter::vdw const & tinker::parameter::parameters::find_vdw (std::size_t type) const
{
  std::size_t const N(m_vdws.size());
  if (indextype(VDW) == GROUP) type = group_by_type(type);
  if (type < N)
  { // we look if type is the correct index and we look at +1 and -1
    if (m_vdws[type].index == type) return m_vdws[type];
    else if (type > 0 && (m_vdws[type-1].index == type)) return m_vdws[type-1];
    else if (type < N-1 && (m_vdws[type+1].index == type)) return m_vdws[type+1];
  }
  for (auto const & i : m_vdws) if (i.index == type) return i;
  throw std::logic_error("No suitable vdw parameter found.");
}

tinker::parameter::vdw const & tinker::parameter::parameters::find_vdw14 (std::size_t itype) const
{
  std::size_t const type((indextype(VDW14) == GROUP) ? group_by_type(itype) : itype);
  for (auto const & vdw : m_vdw14s) 
  {
    if (vdw.index == type) return vdw;
  }
  return find_vdw(itype);
}

tinker::parameter::charge const & tinker::parameter::parameters::find_chg (std::size_t type) const
{
  std::size_t const N(m_charges.size());
  if (indextype(CHARGE) == GROUP) type = group_by_type(type);
  if (type < N)
  { // we look if type is the correct index and we look at +1 and -1
    if (m_charges[type].index == type) return m_charges[type];
    else if (type > 0 && (m_charges[type-1].index == type)) return m_charges[type-1];
    else if (type < N-1 && (m_charges[type+1].index == type)) return m_charges[type+1];
  }
  for (auto const & i : m_charges) if (i.index == type) return i;
  throw std::logic_error("No suitable charge parameter found.");
}


/*

   ###    ########  #######  ##     ## 
  ## ##      ##    ##     ## ###   ### 
 ##   ##     ##    ##     ## #### #### 
##     ##    ##    ##     ## ## ### ## 
#########    ##    ##     ## ##     ## 
##     ##    ##    ##     ## ##     ## 
##     ##    ##     #######  ##     ## 

*/

tinker::parameter::atom::atom (std::string const &line)
  : mass(), type(), group(), atomic(), bonds()
{
  std::size_t const sd(line.find_first_of("\"")+1), lq(line.find_last_of("\"")), ld(lq-sd);
  description = line.substr(sd, ld);
  std::istringstream firstpart(line.substr(4, sd-1));
  firstpart >> type >> group >> symbol;
  std::istringstream secondpart(line.substr(lq+1, line.length()-lq-1));
  secondpart >> atomic >> mass >> bonds;
}

// Calculates the amount of memory (in byte) neccessary
// to store a specific tinker::parameter::atom
std::size_t tinker::parameter::atom::req_mem (void) const
{
  std::size_t symb_size = symbol.capacity()*sizeof(std::string::traits_type::char_type);
  std::size_t desc_size = description.capacity()*sizeof(std::string::traits_type::char_type);
  return sizeof(atom)+symb_size+desc_size;
}

/*

##     ## ########  ######## ##    ## 
##     ## ##     ## ##        ##  ##  
##     ## ##     ## ##         ####   
##     ## ########  ######      ##    
##     ## ##   ##   ##          ##    
##     ## ##    ##  ##          ##    
 #######  ##     ## ########    ##    

*/

tinker::parameter::ureybrad::ureybrad (std::string const &line)
   : f(), ideal(), index()
{
  std::size_t const s(line.find_last_of("d")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[3];
  lss >> indices[0] >> indices[1] >> indices[2] >> f >> ideal;
  index[0] = std::abs(indices[0]);
  index[1] = std::abs(indices[1]);
  index[2] = std::abs(indices[2]);
}

bool tinker::parameter::ureybrad::check (std::size_t a, std::size_t b, std::size_t c) const
{
  if (b == index[1] && ((a == index[0] && c == index[2]) || 
                        (a == index[2] && c == index[0]))) return true;
  return false;
}

/*

   ###    ##    ##  ######   ##       ######## 
  ## ##   ###   ## ##    ##  ##       ##       
 ##   ##  ####  ## ##        ##       ##       
##     ## ## ## ## ##   #### ##       ######   
######### ##  #### ##    ##  ##       ##       
##     ## ##   ### ##    ##  ##       ##       
##     ## ##    ##  ######   ######## ######## 

*/

tinker::parameter::angle::angle (std::string const &line)
   : f(), ideal(), index()
{
  std::size_t const s(line.find_last_of("e")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[3];
  lss >> indices[0] >> indices[1] >> indices[2] >> f >> ideal;
  index[0] = std::abs(indices[0]);
  index[1] = std::abs(indices[1]);
  index[2] = std::abs(indices[2]);
}

bool tinker::parameter::angle::check (std::size_t a, std::size_t b, std::size_t c) const
{
  if (b == index[1] && ((a == index[0] && c == index[2]) || 
                        (a == index[2] && c == index[0]))) return true;
  return false;
}

/*

########   #######  ##    ## ########  
##     ## ##     ## ###   ## ##     ## 
##     ## ##     ## ####  ## ##     ## 
########  ##     ## ## ## ## ##     ## 
##     ## ##     ## ##  #### ##     ## 
##     ## ##     ## ##   ### ##     ## 
########   #######  ##    ## ########  

*/

tinker::parameter::bond::bond (std::string const &line)
   : f(), ideal(), index()
{
  std::size_t const s(line.find_last_of("d")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[2];
  lss >> indices[0] >> indices[1] >> f >> ideal;
  index[0] = std::abs(indices[0]);
  index[1] = std::abs(indices[1]);
}

bool tinker::parameter::bond::check (std::size_t a, std::size_t b) const
{
  if ((a == index[0] && b == index[1]) || (a == index[1] && b == index[0])) return true;
  return false;
}

/*

 ######  ##     ##    ###    ########   ######   ######## 
##    ## ##     ##   ## ##   ##     ## ##    ##  ##       
##       ##     ##  ##   ##  ##     ## ##        ##       
##       ######### ##     ## ########  ##   #### ######   
##       ##     ## ######### ##   ##   ##    ##  ##       
##    ## ##     ## ##     ## ##    ##  ##    ##  ##       
 ######  ##     ## ##     ## ##     ##  ######   ######## 

*/

tinker::parameter::charge::charge (std::string const &line)
   : c(), index()
{
  std::size_t const s(line.find_last_of("e")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[1];
  lss >> indices[0] >> c;
  index = std::abs(indices[0]);
}

bool tinker::parameter::charge::check (std::size_t a) const
{
  return (a == index);
}

/*

########  #######  ########   ######  ####  #######  ##    ## 
   ##    ##     ## ##     ## ##    ##  ##  ##     ## ###   ## 
   ##    ##     ## ##     ## ##        ##  ##     ## ####  ## 
   ##    ##     ## ########   ######   ##  ##     ## ## ## ## 
   ##    ##     ## ##   ##         ##  ##  ##     ## ##  #### 
   ##    ##     ## ##    ##  ##    ##  ##  ##     ## ##   ### 
   ##     #######  ##     ##  ######  ####  #######  ##    ## 

*/

tinker::parameter::torsion::torsion (std::string const &line)
  : force(), ideal(), index(), order(), number(), max_order()
{
  using std::abs;
  std::size_t const s(line.find_first_of(" ")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[4];
  lss >> indices[0] >> indices[1] >> indices[2] >> indices[3];
  index[0] = std::abs(indices[0]);
  index[1] = std::abs(indices[1]);
  index[2] = std::abs(indices[2]);
  index[3] = std::abs(indices[3]);
  for (std::size_t i(0u); i<4; ++i)
  {
    lss >> force[number] >> ideal[number] >> order[number];
    if (abs(force[number]) > 0.0)
    {
      max_order = std::max(max_order, order[number]);
      ++number;
    }
  }
}

bool tinker::parameter::torsion::empty (void) const
{
  return !((fabs(force[0]) > 0.0) || (fabs(force[1]) > 0.0) || (fabs(force[2]) > 0.0) || (fabs(force[3]) > 0.0));
}


std::size_t tinker::parameter::torsion::check2 (std::size_t a, std::size_t b, std::size_t c, std::size_t d) const
{
  std::size_t fit(0u);
  bool valid[4] = {false, false, false, false};
  if(index[0] == a) 
  {
    valid[0] = true;
    ++fit;
  } 
  else if (index[0] == 0) valid[0] = true;
  if(index[1] == b) 
  {
    valid[1] = true;
    ++fit;
  } 
  else if (index[1] == 0) valid[1] = true;
  if(index[2] == c) 
  {
    valid[2] = true;
    ++fit;
  } 
  else if (index[2] == 0) valid[2] = true;
  if(index[3] == d) 
  {
    valid[3] = true;
    ++fit;
  } 
  else if (index[3] == 0) valid[3] = true;
  if (valid[0] && valid[1] && valid[2] && valid[3]) return fit;
  return 0u;
}

/**looks for torsion
returns 4 if parameters in parameter file only with real atoms
returns smaller numbers if some of the atom types in paramater file are 0
returns 0 if not fitting in*/
std::size_t tinker::parameter::torsion::check (std::size_t a, std::size_t b, std::size_t c, std::size_t d) const
{

  if (index[0] == a && index[1] == b && index[2] == c && index[3] == d) return 4u;
  else if (index[0] == a && index[1] == b && index[2] == c && index[3] == d) return 4u;
  else return std::max(check2(a,b,c,d), check2(d,c,b,a));

}

/*
 Imptors
####         ########  #######  ########   ######  ####  #######  ##    ## 
 ##             ##    ##     ## ##     ## ##    ##  ##  ##     ## ###   ## 
 ##             ##    ##     ## ##     ## ##        ##  ##     ## ####  ## 
 ##  #######    ##    ##     ## ########   ######   ##  ##     ## ## ## ## 
 ##             ##    ##     ## ##   ##         ##  ##  ##     ## ##  #### 
 ##             ##    ##     ## ##    ##  ##    ##  ##  ##     ## ##   ### 
####            ##     #######  ##     ##  ######  ####  #######  ##    ## 

 const std::size_t tinker::parameter::global::improper_sequence[4] = 
  {tinker::parameter::global::CENTER, 
   tinker::parameter::global::LIG1, 
   tinker::parameter::global::LIG2, 
   tinker::parameter::global::TWIST};

const std::size_t tinker::parameter::global::imptors_sequence[4] = 
  {tinker::parameter::global::LIG1, 
   tinker::parameter::global::LIG2, 
   tinker::parameter::global::CENTER, 
   tinker::parameter::global::TWIST};

*/

tinker::parameter::imptor::imptor (std::string const &line)
  : force(), ideal(), order(), number(), max_order(), center(), ligand()
{
  using std::abs;
  std::size_t const s(line.find_first_of(" ")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[4];
  lss >> indices[0] >> indices[1] >> indices[2] >> indices[3];
  center = std::abs(indices[2]);
  ligand[0] = std::abs(indices[0]);
  ligand[1] = std::abs(indices[1]);
  ligand[2] = std::abs(indices[3]);
  for (std::size_t i(0u); i<4; ++i)
  {
    lss >> force[number] >> ideal[number] >> order[number];
    if (abs(force[number]) > 0.0)
    {
      max_order = std::max(max_order, order[number]);
      ++number;
    }
  }
}

bool tinker::parameter::imptor::empty (void) const
{
  using std::abs;
  return !((abs(force[0]) > 0.0) || (abs(force[1]) > 0.0) || (abs(force[2]) > 0.0));
}

/*
 Improper
####         ########  #######  ########   ######  ####  #######  ##    ## 
 ##             ##    ##     ## ##     ## ##    ##  ##  ##     ## ###   ## 
 ##             ##    ##     ## ##     ## ##        ##  ##     ## ####  ## 
 ##  #######    ##    ##     ## ########   ######   ##  ##     ## ## ## ## 
 ##             ##    ##     ## ##   ##         ##  ##  ##     ## ##  #### 
 ##             ##    ##     ## ##    ##  ##    ##  ##  ##     ## ##   ### 
####            ##     #######  ##     ##  ######  ####  #######  ##    ## 

*/

tinker::parameter::improper::improper (std::string const &line)
  : force(), ideal(), order(), number(), max_order(), center(), ligand()
{
  std::size_t const s(line.find_first_of(" ")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[4];
  lss >> indices[0] >> indices[1] >> indices[2] >> indices[3];
  center = std::abs(indices[0]);
  ligand[0] = std::abs(indices[3]);
  ligand[1] = std::abs(indices[1]);
  ligand[2] = std::abs(indices[2]);
  for (std::size_t i(0u); i<4; ++i)
  {
    lss >> force[number] >> ideal[number] >> order[number];
    if (fabs(force[number]) > 0.0)
    {
      max_order = std::max(max_order, order[number]);
      ++number;
    }
  }
}

bool tinker::parameter::improper::empty (void) const
{
  return !((fabs(force[0]) > 0.0) || (fabs(force[1]) > 0.0) || (fabs(force[2]) > 0.0));
}

/*

##     ## ##     ## ##       ######## #### ########   #######  ##       ######## 
###   ### ##     ## ##          ##     ##  ##     ## ##     ## ##       ##       
#### #### ##     ## ##          ##     ##  ##     ## ##     ## ##       ##       
## ### ## ##     ## ##          ##     ##  ########  ##     ## ##       ######   
##     ## ##     ## ##          ##     ##  ##        ##     ## ##       ##       
##     ## ##     ## ##          ##     ##  ##        ##     ## ##       ##       
##     ##  #######  ########    ##    #### ##         #######  ######## ########

*/

tinker::parameter::multipole::multipole(std::string const &line)
: charge(), dipole(), quadrupole(), index(), has_4th_index()
{
	std::size_t const s(line.find_first_of("e") + 1), l(line.length() - s);
	std::istringstream lss(line.substr(s, l));
  ptrdiff_t indices[4] = { 0 };
	//double d2(0.0);
	//ptrdiff_t d1(0);
	axt = Z_THEN_X;
  std::vector<double> tmp;
  double x;
  while (lss >> x)
  {
    tmp.push_back(x);
  }
  if (tmp.size() == 4u)
  {
    charge = tmp[3u];
    indices[0u] = static_cast<ptrdiff_t>(tmp[0u]);
    indices[1u] = static_cast<ptrdiff_t>(tmp[1u]);
    indices[2u] = static_cast<ptrdiff_t>(tmp[2u]);
  }
  if (tmp.size() == 5u)
  {
    charge = tmp[4u];
    indices[0u] = static_cast<ptrdiff_t>(tmp[0u]);
    indices[1u] = static_cast<ptrdiff_t>(tmp[1u]);
    indices[2u] = static_cast<ptrdiff_t>(tmp[2u]);
    indices[3u] = static_cast<ptrdiff_t>(tmp[3u]);
  }
	/*lss >> indices[0] >> indices[1] >> indices[2];
  if (!(lss >> indices[3]) || !(lss >> charge))
  {
    lss.str(line.substr(s, l));
    lss >> indices[0] >> indices[1] >> indices[2] >> charge;
  }*/
  /*
	sscanf(lss.str().c_str(), " %d %d %d %d %lf", &indices[0], &indices[1], &indices[2], &indices[3], &charge);
	if (indices[3] == 0){
		sscanf(lss.str().c_str(), "%d %d %d %lf", &indices[0], &indices[1], &indices[2], &charge);
	}*/
	index[0] = abs(indices[0]);
	index[1] = abs(indices[1]);
	index[2] = abs(indices[2]);
	//indices[3] = d1;
	index[3] = abs(indices[3]);
	//charge = d2;
	if (indices[1] == 0) axt = NONE;
	else if (indices[1] != 0 && indices[2] == 0) axt = Z_AXIS;
	else if (indices[1] <= 0 || indices[2] <= 0) axt = BISECTOR;
	else if (indices[2] <= 0 && indices[3] <= 0) axt = Z_BISECTOR;
	else if (std::max(std::max(indices[1], indices[2]), indices[3]) < 0) axt = THREEFOLD;
}

void tinker::parameter::multipole::parse_line_2(std::string const & line)
{
	std::istringstream lss(line);
	lss >> dipole.x() >> dipole.y() >> dipole.z();
}
void tinker::parameter::multipole::parse_line_3(std::string const & line)
{
	std::istringstream lss(line);
	lss >> quadrupole[0u];
}
void tinker::parameter::multipole::parse_line_4(std::string const & line)
{
	std::istringstream lss(line);
	lss >> quadrupole[1u] >> quadrupole[2u];
}
void tinker::parameter::multipole::parse_line_5(std::string const & line)
{
	std::istringstream lss(line);
	lss >> quadrupole[3u] >> quadrupole[4u] >> quadrupole[5u];
}
void tinker::parameter::multipole::multiply (double const quadro_factor, double const dip_factor)
{
	dipole *= dip_factor;
	quadrupole *= quadro_factor;
}

/*

 #######  ########  ########  ######## ##    ## ########  
##     ## ##     ## ##     ## ##       ###   ## ##     ## 
##     ## ##     ## ##     ## ##       ####  ## ##     ## 
##     ## ########  ########  ######   ## ## ## ##     ## 
##     ## ##        ##     ## ##       ##  #### ##     ## 
##     ## ##        ##     ## ##       ##   ### ##     ## 
 #######  ##        ########  ######## ##    ## ########    

*/

tinker::parameter::opbend::opbend (std::string const &line)
  : f(), index()
{
  std::size_t const s(line.find_first_of("d")+1), l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t indices[2];
  double d1(0.0), d2(0.0);
  lss >> indices[0] >> indices[1] >> f >> d1 >> d2;
  index[0] = std::abs(indices[0]);
  index[1] = std::abs(indices[1]);
  if (std::fabs(d2) > 0.0)
  {
    index[2] = static_cast<std::size_t>(std::fabs(f));
    index[3] = static_cast<std::size_t>(std::fabs(d1));
    std::swap(f,d2);
  }
  else if (std::fabs(d1) > 0.0)
  {
    index[2] = static_cast<std::size_t>(std::fabs(f));
    std::swap(f,d1);
  }
}

bool tinker::parameter::opbend::check(std::size_t a, std::size_t b) const
{

	if ((index[0] == a && index[1] == b) || (index[0] == b && index[1] == a)) return true;
	return false;

}
/*

########   #######  ##          ###    ########  #### ######## ######## 
##     ## ##     ## ##         ## ##   ##     ##  ##       ##  ##       
##     ## ##     ## ##        ##   ##  ##     ##  ##      ##   ##       
########  ##     ## ##       ##     ## ########   ##     ##    ######   
##        ##     ## ##       ######### ##   ##    ##    ##     ##       
##        ##     ## ##       ##     ## ##    ##   ##   ##      ##       
##         #######  ######## ##     ## ##     ## #### ######## ######## 

*/

tinker::parameter::polarize::polarize(std::string const &line)
: p(), pp(), index(),bonded()
{
	std::size_t const s(line.find_first_of("e") + 1), l(line.length() - s);
	std::size_t const n_fl = std::count(line.begin(), line.end(), '.');
	std::istringstream lss(line.substr(s, l));
	ptrdiff_t indices[4] = { 0, 0, 0, 0 };
	//bond_index.resize(3U);
	lss >> indices[0];
	if (n_fl == 2)  lss >> p >> pp;
	else lss >> p;

	index = std::abs(indices[0]);

	ptrdiff_t x;
	std::vector<ptrdiff_t> tmp;
	while (lss >> x)
	{
		tmp.push_back(x);
	}
	if (tmp.size() == 1u)
	{
		
		indices[1] = std::abs(tmp[0U]);
	
	}
	if (tmp.size() == 2u)
	{
		indices[1] = std::abs(tmp[0U]);
		indices[2] = std::abs(tmp[1U]);
		
	}
	if (tmp.size() == 3u)
	{
		indices[1] = std::abs(tmp[0U]);
		indices[2] = std::abs(tmp[1U]);
		indices[3] = std::abs(tmp[2U]);

	}
	bonded[0] = std::abs(indices[1]);
	bonded[1] = std::abs(indices[2]);
	bonded[2] = std::abs(indices[3]);
	not_contracted_bond[0] = bonded[0];
	not_contracted_bond[1] = bonded[1];
	not_contracted_bond[2] = bonded[2];

	
}

/*

 ######  ######## ########  ########  ######## ##    ## ########  
##    ##    ##    ##     ## ##     ## ##       ###   ## ##     ## 
##          ##    ##     ## ##     ## ##       ####  ## ##     ## 
 ######     ##    ########  ########  ######   ## ## ## ##     ## 
      ##    ##    ##   ##   ##     ## ##       ##  #### ##     ## 
##    ##    ##    ##    ##  ##     ## ##       ##   ### ##     ## 
 ######     ##    ##     ## ########  ######## ##    ## ########  

*/

tinker::parameter::strbend::strbend (std::string const &line)
  : f(), index()
{
  std::size_t const s(line.find_first_of("d")+1), l(line.length()-s);
  std::string const sub(line.substr(s, l));
  //std::size_t const num_floats = std::count(sub.begin(), sub.end(), '.');
  std::istringstream lss(sub);
  /*double d[6] = {};*/
  lss >> index[0] >> index[1] >> index[2] >> f >> ff;
  //std::size_t last_valid(0);
  //for (std::size_t j=0; j<6; ++j) 
  //{
  //  last_valid = 5-j;
  //  if (fabs(d[5-j]) > 0.0) break;
  //}
  //std::size_t const num_indices(last_valid+1-num_floats);
  //for (std::size_t i(0); i<=last_valid; ++i)
  //{
  //  if (i<num_indices) index[i] = static_cast<std::size_t>(std::fabs(d[i]));
  //  else f[i-num_indices] = d[i];
  //}
}


std::size_t tinker::parameter::strbend::check(std::size_t a, std::size_t b, std::size_t c) const
{

	if (b == index[1] && ((a == index[0] && c == index[2]) ||
		(a == index[2] && c == index[0])))
	{
		if (a <= c)return 0U;
		else return 1U;
	}
	else return 2U;
}

/*

##     ## ########  ##      ## 
##     ## ##     ## ##  ##  ## 
##     ## ##     ## ##  ##  ## 
##     ## ##     ## ##  ##  ## 
 ##   ##  ##     ## ##  ##  ## 
  ## ##   ##     ## ##  ##  ## 
   ###    ########   ###  ###  

*/

tinker::parameter::vdw::vdw (std::string const &line)
   : r(), e(), rf(), index()
{
  std::size_t s(line.find_last_of("w")+1);
  if (line[s] == '1') s+=2;
  std::size_t const l(line.length()-s);
  std::istringstream lss(line.substr(s, l));
  std::ptrdiff_t ind;

  lss >> ind >> r >> e >> rf;
  index = std::abs(ind);
  
 

}

bool tinker::parameter::vdw::check(size_t a, size_t b) const
{

	if ((indices[0] == a && indices[1] == b) || (indices[0] == b && indices[1] == a)) return true;
	return false;

}

bool tinker::parameter::vdw::check (size_t a) const
{
  return (a == index);
}

/*

######## #### ##       ########            ####          #######  
##        ##  ##       ##                   ##          ##     ## 
##        ##  ##       ##                   ##          ##     ## 
######    ##  ##       ######               ##  ####### ##     ## 
##        ##  ##       ##                   ##          ##     ## 
##        ##  ##       ##                   ##          ##     ## 
##       #### ######## ########            ####          #######  
 
*/

/*

 ######   ######## ##    ## ######## ########     ###    ##       
##    ##  ##       ###   ## ##       ##     ##   ## ##   ##       
##        ##       ####  ## ##       ##     ##  ##   ##  ##       
##   #### ######   ## ## ## ######   ########  ##     ## ##       
##    ##  ##       ##  #### ##       ##   ##   ######### ##       
##    ##  ##       ##   ### ##       ##    ##  ##     ## ##       
 ######   ######## ##    ## ######## ##     ## ##     ## ######## 

*/

void tinker::parameter::forcefield_types::parse_line (std::string const & line) 
{
  if      (line.find("OPLS") != std::string::npos) { value = OPLSAA; }
  else if (line.find("CHARMM") != std::string::npos) { value = CHARMM22; }
  else if (line.find("AMBER") != std::string::npos) { value = AMBER; }
  else if (line.find("AMOEBA") != std::string::npos) { value = AMOEBA; } 
  else throw std::runtime_error(std::string("Forcefield is unknown: ").append(line));
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, forcefield_types const & t)
{
  if (t.value==forcefield_types::OPLSAA) stream << "OPLS-AA";
  else if (t.value==forcefield_types::CHARMM22) stream << "CHARMM22";
  else if (t.value==forcefield_types::AMBER) stream << "AMBER-FF99";
  else if (t.value==forcefield_types::AMOEBA) stream << "AMOEBA-2009";
  return stream;
}

void tinker::parameter::vdw_types::parse_line (std::string const & line) 
{
  if      (line.find("LENNARD-JONES") != std::string::npos) { value = LENNARD_JONES; }
  else if (line.find("BUFFERED-14-7") != std::string::npos) { value = BUFFERED_14_7; }
  else if (line.find("BUCKINGHAM") != std::string::npos) { value = BUCKINGHAM; }
  else if (line.find("MM3-HBOND") != std::string::npos) { value = MM3_HBOND; }
  else if (line.find("GAUSSIAN") != std::string::npos) { value = GAUSSIAN; }
  else throw std::runtime_error(std::string("Vdwtype is unknow: n").append(line));
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, vdw_types const & t)
{
  if (t.value==vdw_types::LENNARD_JONES) stream << "LENNDARD-JONES";
  else if (t.value==vdw_types::BUFFERED_14_7) stream << "BUFFERED-14-7";
  else if (t.value==vdw_types::BUCKINGHAM) stream << "BUCKINGHAM";
  else if (t.value==vdw_types::MM3_HBOND) stream << "MM3-HBOND";
  else if (t.value==vdw_types::GAUSSIAN) stream << "GAUSSIAN";
  return stream;
}

void tinker::parameter::average_rules::parse_line (std::string const & line) 
{
  if      (line.find("GEOMETRIC") != std::string::npos) { value = GEOMETRIC; }
  else if (line.find("ARITHMETIC") != std::string::npos) { value = ARITHMETIC; }
  else if (line.find("CUBIC-MEAN") != std::string::npos) { value = CUBIC_MEAN; }
  else if (line.find("HARMONIC") != std::string::npos) { value = HARMONIC; }
  else if (line.find("HHG") != std::string::npos) { value = HHG; }
  else throw std::runtime_error(std::string("Avarage rule is unknown: ").append(line));
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, average_rules const & t)
{
  if (t.value==average_rules::GEOMETRIC) stream << "GEOMETRIC";
  else if (t.value==average_rules::ARITHMETIC) stream << "ARITHMETIC";
  else if (t.value==average_rules::CUBIC_MEAN) stream << "CUBIC-MEAN";
  return stream;
}

void tinker::parameter::radius_types::parse_line (std::string const & line) 
{
  if      (line.find("SIGMA") != std::string::npos) { value = SIGMA; }
  else if (line.find("R-MIN") != std::string::npos) { value = R_MIN; }
  else throw std::runtime_error(std::string("Radius Type is unknown: ").append(line));
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, radius_types const & t)
{
  if (t.value==radius_types::SIGMA) stream << "SIGMA";
  else if (t.value==radius_types::R_MIN) stream << "R-MIN";
  return stream;
}

void tinker::parameter::radius_sizes::parse_line (std::string const & line) 
{
  if      (line.find("DIAMETER") != std::string::npos) { value = DIAMETER; }
  else if (line.find("RADIUS") != std::string::npos) { value = RADIUS; }
  else throw std::runtime_error(std::string("Radius Size is unknown: ").append(line));
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, radius_sizes const & t)
{
  if (t.value==radius_sizes::RADIUS) stream << "RADIUS";
  else if (t.value==radius_sizes::DIAMETER) stream << "DIAMETER";
  return stream;
}

void tinker::parameter::scales::parse_line (std::string const & line) 
{
  std::size_t const s1(line.find_first_of("-")+2), l1(line.find_last_of("-")-s1);
  std::size_t i;
  std::istringstream(line.substr(s1, l1)) >> i; // 1 for 11, 2 for 12, 3 for 13...
  --i; // 0 for 11, 1 for 12, ---
  std::size_t const s2(line.find_last_of("e")+1), l2(line.length()-s2);
  std::istringstream(line.substr(s2, l2)) >> value[i];
  if (fabs(value[i]) > 0.0) use[i] = true;
  else use[i] = false;
}

void tinker::parameter::scales::to_stream (std::ostream & stream, std::string type) const
{
  for (std::size_t i(0u); i<5u; ++i)
  {
    if (use[i])
    {
      std::stringstream ss;
      ss << type << "-1" << i+1 << "-scale";
      stream << std::left << std::setw(30) << ss.str() << value[i];
      stream << std::endl;
    }
  }
}

void tinker::parameter::factors::parse_line (std::string const & line) 
{
  std::size_t i;
  if      (line.find("-cubic") != std::string::npos) 
  { 
    i = 0; 
    order = std::max(order, std::size_t(1u)); 
  }
  else if (line.find("-quartic") != std::string::npos) 
  { 
    i = 1; 
    order = std::max(order, std::size_t(2u)); 
  }
  else if (line.find("-pentic") != std::string::npos) 
  { 
    i = 2; 
    order = std::max(order, std::size_t(3u)); 
  }
  else if (line.find("-sextic") != std::string::npos) 
  { 
    i = 3; 
    order = std::max(order, std::size_t(4u)); 
  }
  else  throw std::runtime_error(std::string("Factor of unknown order: ").append(line));
  std::istringstream(line.substr(line.find_last_of("c")+1, line.length())) >> value[i];
}

void tinker::parameter::index::parse_line (std::string const & line) 
{
  if (line.find("TYPE") != std::string::npos) value = TYPE;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, index const & t)
{
  if (t.value==index_types::GROUP) stream << "CLASS";
  else if (t.value==index_types::TYPE) stream << "TYPE";
  return stream;
}

/**
 * This functions reads the ALL forcefield parameters
 * from a param-file.
 *
 * @param filename: filename of the parameter file (needs to be in the same folder as executable)
 */
void tinker::parameter::parameters::from_file (std::string const & filename)
{
  // LBL_FileReader reads a file and puts content into the public member "data"
  LBL_FileReader param_file(filename);

  // If data is empty, throw error
  if (param_file.data.empty())
  {
    std::string fn("'");
    fn.append(filename);
    fn.append("'");
    throw std::runtime_error(std::string("File for TINKER-like parameters not found at ").append(fn));
  }

  // Store all the TINKER parameters from the para_file
  // In the parameters object ("this")
  parse_lines(param_file.data);


  if (!m_atoms.empty())
  {
    std::sort(m_atoms.begin(), m_atoms.end(), atom::typeless());
    std::size_t const maxtype(m_atoms.back().type);

    std::vector<atom> m_atoms_temporary(maxtype);

    // At first I thought: Appareantly we adjust for zero or one based indexation
    // However, after closer examination, we actually do NOTHING o.O
    // We create an equal vector by obscure means
    for (auto i : m_atoms) m_atoms_temporary[i.type - 1] = i;
    // In my current setup, they are euqal by the test:
    //bool testIfEqual = true;
    //for (size_t i = 0u; i < maxtype; i++)
    //if (!(m_atoms[i] == m_atoms_temporary[i]))
    //{
    //  testIfEqual = false;
    //}
    // However, maybe there is a reason why this is here.
    // Maybe it is necessary for other forcefields
    // I dont know...


    m_atoms = m_atoms_temporary;
  } 
  m_valid = true;
}

// + HELPERS ##########
static inline bool str_match_begin (std::string const &line, std::string const &check)
{
  return (line.substr(0, check.length()) == check);
}

template<class T>
static inline bool check_string_put (std::string const &line, std::string const &check, T & target)
{
  if (str_match_begin(line, check)) 
  {
    std::istringstream(line.substr(check.length(),line.length()-check.length())) >> target;
    return true;
  } else 
    return false;
}

/*

######## #### ##       ########            ####          #######  
##        ##  ##       ##                   ##          ##     ## 
##        ##  ##       ##                   ##          ##     ## 
######    ##  ##       ######               ##  ####### ##     ## 
##        ##  ##       ##                   ##          ##     ## 
##        ##  ##       ##                   ##          ##     ## 
##       #### ######## ########            ####          #######  
 
*/

/*

########     ###    ########   ######  ######## 
##     ##   ## ##   ##     ## ##    ## ##       
##     ##  ##   ##  ##     ## ##       ##       
########  ##     ## ########   ######  ######   
##        ######### ##   ##         ## ##       
##        ##     ## ##    ##  ##    ## ##       
##        ##     ## ##     ##  ######  ######## 

*/

/**
 * This functions reads the forcefield parameters the lines of a param-file.
 * These "lines" are obtained using LBL_FileReader from filemanipulation.h
 * The data is then properly stored in the members of the parameters object.
 * Yes, ALL the parameters from the file are stored, even the parameters of
 * atoms that are maybe not needed.
 *
 * @param lines: vector containing the lines of a properly formated file. Usually obtained using
 * LBL_FileReader object from filemanipulation.h
 */
void tinker::parameter::parameters::parse_lines (std::vector<std::string> const & lines)
{
  std::size_t const numberOfLines(lines.size());
  for (std::size_t i(0u); i <  numberOfLines; ++i)
  {
    std::string const & line(lines[i]);
    if (line.empty()) continue;
    switch (line[0u])
    {
    case 'a':
      { // angle, atom
        if (str_match_begin(line, "atom "))  m_atoms.push_back(line);
        else if (str_match_begin(line, "angle ")) m_angles.push_back(line);
        else if (str_match_begin(line, "angle-")) m_general.angle_factor.parse_line(line);
        else check_string_put(line, "angleunit", m_general.angleunit);
        break;
      }
    case 'b':
      {
        if (str_match_begin(line, "bond ")) m_bonds.push_back(line);
        else if(str_match_begin(line, "bond-")) m_general.bond_factor.parse_line(line);
        else check_string_put(line, "bondunit", m_general.bondunit);
        break;
      }
    case 'c':
      {
        if (str_match_begin(line, "charge ")) m_charges.push_back(line);
        else if (str_match_begin(line, "chg-")) m_general.chg_scale.parse_line(line);
        break;
      }
    case 'd':
      {
        if (!check_string_put(line, "dielectric", m_general.dielectric))
        {
          if (str_match_begin(line, "direct-")) m_general.direct_scale.parse_line(line);
        }
        break;
      }
    case 'e':
      {
        if (!check_string_put(line, "electric", m_general.electric))
        {
          if (str_match_begin(line, "epsilonrule")) m_general.epsilonrule.parse_line(line);
        }
        break;
      }
    case 'f':
      {
        if (str_match_begin(line, "forcefield")) m_general.forcefield.parse_line(line);
        break;
      }
    case 'i':
      {
        if (!check_string_put(line, "imptorunit", m_general.imptorunit))
        {
          if (str_match_begin(line, "imptors")) m_imptors.push_back(line);
          else if (str_match_begin(line, "improper")) m_impropers.push_back(line);
        }
        break;
      }
    case 'm':
      {
        if (str_match_begin(line, "multipole")) 
        {
          m_multipoles.push_back(line);
          m_multipoles.back().parse_line_2(lines[++i]);
          m_multipoles.back().parse_line_3(lines[++i]);
          m_multipoles.back().parse_line_4(lines[++i]);
          m_multipoles.back().parse_line_5(lines[++i]);
        }
        else if (str_match_begin(line, "mpole-")) m_general.mpole_scale.parse_line(line);
        else if (str_match_begin(line, "mutual-")) m_general.mutual_scale.parse_line(line);
        break;
      }
    case 'o':
      {
		if (str_match_begin(line, "opbend ")) m_opbends.push_back(line); 
        else if (str_match_begin(line, "opbend-")) m_general.opbend_factor.parse_line(line);
        break;
      }
    case 'p':
      {
        if (str_match_begin(line, "polarize ")) m_polarizes.push_back(line);
        else if (str_match_begin(line, "polar-")) m_general.polar_scale.parse_line(line);
        break;
      }
    case 'r':
      {
        if (str_match_begin(line, "radiustype")) m_general.radiustype.parse_line(line);
        else if (str_match_begin(line, "radiussize")) m_general.radiussize.parse_line(line);
        else if (str_match_begin(line, "radiusrule")) m_general.radiusrule.parse_line(line);
        
        break;
      }
    case 's':
      {
				if (str_match_begin(line, "strbnd ")) m_strbends.push_back(line);
        break;
      }
    case 't':
      {
        if (!check_string_put(line, "torsionunit", m_general.torsionunit))
        {
          if (str_match_begin(line, "torsion ")) m_torsions.push_back(line);
        }
        break;
      }
    case 'u':
      {
        if (str_match_begin(line, "ureybrad ")) m_ureybrads.push_back(line);
        break;
      }
    case 'v':
      {
        if (str_match_begin(line, "vdw14 ")) m_vdw14s.push_back(line);
        else if (str_match_begin(line, "vdw-")) m_general.vdw_scale.parse_line(line);
        else if (str_match_begin(line, "vdwtype")) m_general.vdwtype.parse_line(line);
        else if (str_match_begin(line, "vdwindex")) m_general.indices[VDW].parse_line(line);
        else if (str_match_begin(line, "vdw ")) m_vdws.push_back(line);
        break;
      }
    }
  }
}

/*

 ######   #######  ##    ## ######## ########     ###     ######  ######## ####  #######  ##    ## 
##    ## ##     ## ###   ##    ##    ##     ##   ## ##   ##    ##    ##     ##  ##     ## ###   ## 
##       ##     ## ####  ##    ##    ##     ##  ##   ##  ##          ##     ##  ##     ## ####  ## 
##       ##     ## ## ## ##    ##    ########  ##     ## ##          ##     ##  ##     ## ## ## ## 
##       ##     ## ##  ####    ##    ##   ##   ######### ##          ##     ##  ##     ## ##  #### 
##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##    ##    ##     ##  ##     ## ##   ### 
 ######   #######  ##    ##    ##    ##     ## ##     ##  ######     ##    ####  #######  ##    ## 

*/


/*! Contract the TINKER atom types found in the .arc TINKER coordinate inputfile.
 */
tinker::parameter::parameters tinker::parameter::parameters::contract(std::vector<size_t> actual_types) const
{
  std::sort(actual_types.begin(), actual_types.end(), std::less<std::size_t>());
  // NA: Number of atoms in structure
  // ACT: Number of (actual) unqiue forcefield atom types plus one
  std::size_t const numberOfAtoms(m_atoms.size()), ACT(actual_types.size()+1);

  // During CASTs initializing procedure an empty coords object is created
  // During this procedure, an empty energy interface is created according to inputfile specifications
  // This following if statement takes care of this for the initial startup procedure.
  // Then, for now, empty parameter container is returned.
  parameters parametersToBeContracted;
  if (actual_types.empty()) return parametersToBeContracted;

  // Tell the parameter container that it describes contracted parameters.
  // It is (appareantly) also possible to work with not-contracted parameters
  parametersToBeContracted.contracted = true;

  // Copy the "global" m_general struct containing assorted stuff
  // I am not sure what the stuff that is in there does
  // But, oh well, guess we copy it to the parameter file
  parametersToBeContracted.m_general = m_general;

  // OK, so our parameter contains private members m_atoms 
  // and m_uncontracted_atoms. Seems reasonable.
  // Adjusting the size of m_atoms, probably actually standing for
  // "Contracted atom types"
  // I think m_atoms are still uncontracted atom types?
  parametersToBeContracted.m_atoms.resize(actual_types.size());
  parametersToBeContracted.m_uncontracted_atoms = m_atoms;
  // Huh, didnt expect this. Interesting

  // Ok, i removed this line. No f*cking clue what its supposed to do
  //enum { GROUP=GROUP, TYPE=TYPE};

  // FIll with "numberOfAtoms" copies of the value of "ACT+1"
  parametersToBeContracted.contraction_map.assign(numberOfAtoms +1, ACT+1);
  parametersToBeContracted.group_contraction.assign(numberOfAtoms +1, ACT+1);

  parametersToBeContracted.contraction_map[0u] = 0u;
  parametersToBeContracted.group_contraction[0u] = 0u;
  parametersToBeContracted.m_reduced_groups.push_back(0u);
  parametersToBeContracted.m_reduced_types.push_back(0u);

  // contract types & groups
  for (auto i : actual_types)
  {
    if (i > m_atoms.size() || i <= 0)
    {
      std::stringstream tmpss;
      tmpss << i << " is not a valid type in parameter file.";
      throw std::runtime_error(tmpss.str().c_str());
    }
    scon::sorted::insert_unique(parametersToBeContracted.m_reduced_groups, m_atoms[i-1].group);
    parametersToBeContracted.contraction_map[i] = parametersToBeContracted.m_reduced_types.size();
    parametersToBeContracted.m_reduced_types.push_back(m_atoms[i-1].type);
  }
  // map groups
  std::size_t const NR(parametersToBeContracted.m_reduced_groups.size());
  for (std::size_t i(0u); i<NR; ++i)
  {
    parametersToBeContracted.group_contraction[parametersToBeContracted.m_reduced_groups[i]] = i;
  }
  // map atom types
  for (auto i : actual_types)
  {
    std::size_t const red_type(parametersToBeContracted.contraction_map[i]),
      red_group(parametersToBeContracted.group_contraction[m_atoms[i-1].group]);
    parametersToBeContracted.m_atoms[red_type-1] = m_atoms[i-1];
    parametersToBeContracted.m_atoms[red_type-1].group = red_group;
    parametersToBeContracted.m_atoms[red_type-1].type = red_type;
  }
  // transfer vdw and charges, 
  for (auto const &i : m_vdws)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index, VDW));
    if (a < ACT)
    {
      parametersToBeContracted.m_vdws.push_back(i);
      parametersToBeContracted.m_vdws.back().index = a;
    }
  }

  for (auto const &i : m_vdw14s)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index, VDW14));
    if (a < ACT)
    {
      parametersToBeContracted.m_vdw14s.push_back(i);
      parametersToBeContracted.m_vdw14s.back().index = a;
    }
  }

  for (auto const &i : m_charges)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index, CHARGE));
    if (a < ACT)
    {
      parametersToBeContracted.m_charges.push_back(i);
      parametersToBeContracted.m_charges.back().index = a;
    }
  }

  for (auto const &i : m_angles)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], ANGLE)),
      b(parametersToBeContracted.contracted_type(i.index[1], ANGLE)),
      c(parametersToBeContracted.contracted_type(i.index[2], ANGLE));
    if ((a < ACT) && (b < ACT) && (c < ACT))
    {
      parametersToBeContracted.m_angles.push_back(i);
      parametersToBeContracted.m_angles.back().index[0] = a;
      parametersToBeContracted.m_angles.back().index[1] = b;
      parametersToBeContracted.m_angles.back().index[2] = c;
    }
  }

  for (auto const &i : m_bonds)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], BOND)),
      b(parametersToBeContracted.contracted_type(i.index[1], BOND));
    if ((a < ACT) && (b < ACT))
    {
      parametersToBeContracted.m_bonds.push_back(i);
      parametersToBeContracted.m_bonds.back().index[0] = a;
      parametersToBeContracted.m_bonds.back().index[1] = b;
    }
  }

  for (auto const &i : m_impropers)
  {
    if (i.empty()) continue;
    std::size_t const a(parametersToBeContracted.contracted_type(i.center, IMPROPER)),
      b(parametersToBeContracted.contracted_type(i.ligand[0], IMPROPER)),
      c(parametersToBeContracted.contracted_type(i.ligand[1], IMPROPER)),
      d(parametersToBeContracted.contracted_type(i.ligand[2], IMPROPER));
    if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT))
    {
      parametersToBeContracted.m_impropers.push_back(i);
      parametersToBeContracted.m_impropers.back().center = a;
      parametersToBeContracted.m_impropers.back().ligand[0] = b;
      parametersToBeContracted.m_impropers.back().ligand[1] = c;
      parametersToBeContracted.m_impropers.back().ligand[2] = d;
    }
  }

  for (auto const &i : m_imptors)
  {
    if (i.empty()) continue;
    std::size_t const a(parametersToBeContracted.contracted_type(i.center, IMPROPER)),
      b(parametersToBeContracted.contracted_type(i.ligand[0], IMPROPER)),
      c(parametersToBeContracted.contracted_type(i.ligand[1], IMPROPER)),
      d(parametersToBeContracted.contracted_type(i.ligand[2], IMPROPER));
    if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT))
    {
      parametersToBeContracted.m_imptors.push_back(i);
      parametersToBeContracted.m_imptors.back().center = a;
      parametersToBeContracted.m_imptors.back().ligand[0] = b;
      parametersToBeContracted.m_imptors.back().ligand[1] = c;
      parametersToBeContracted.m_imptors.back().ligand[2] = d;
    }
  }

  for (auto const &i : m_torsions)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], TORSION)),
      b(parametersToBeContracted.contracted_type(i.index[1], TORSION)),
      c(parametersToBeContracted.contracted_type(i.index[2], TORSION)),
      d(parametersToBeContracted.contracted_type(i.index[3], TORSION));
    if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT))
    {
      parametersToBeContracted.m_torsions.push_back(i);
      parametersToBeContracted.m_torsions.back().index[0] = a;
      parametersToBeContracted.m_torsions.back().index[1] = b;
      parametersToBeContracted.m_torsions.back().index[2] = c;
      parametersToBeContracted.m_torsions.back().index[3] = d;
    }
  }

  for (auto const &i : m_multipoles)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], MULTIPOLE)),
      b(parametersToBeContracted.contracted_type(i.index[1], MULTIPOLE)),
      c(parametersToBeContracted.contracted_type(i.index[2], MULTIPOLE)),
      d(parametersToBeContracted.contracted_type(i.index[3], MULTIPOLE));
    if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT))
    {
      parametersToBeContracted.m_multipoles.push_back(i);
      parametersToBeContracted.m_multipoles.back().index[0] = a;
      parametersToBeContracted.m_multipoles.back().index[1] = b;
      parametersToBeContracted.m_multipoles.back().index[2] = c;
      parametersToBeContracted.m_multipoles.back().index[3] = d;
    }
  }

  for (auto const &i : m_opbends)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], OPBEND)),
      b(parametersToBeContracted.contracted_type(i.index[1], OPBEND)),
      c(parametersToBeContracted.contracted_type(i.index[2], OPBEND)),
      d(parametersToBeContracted.contracted_type(i.index[3], OPBEND));
    if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT))
    {
      parametersToBeContracted.m_opbends.push_back(i);
      parametersToBeContracted.m_opbends.back().index[0] = a;
      parametersToBeContracted.m_opbends.back().index[1] = b;
      parametersToBeContracted.m_opbends.back().index[2] = c;
      parametersToBeContracted.m_opbends.back().index[3] = d;
    }
  }

  for (auto const &i : m_polarizes)
  {
	  size_t const a(parametersToBeContracted.contracted_type(i.index, POLARIZE)),
		  b(parametersToBeContracted.contracted_type(i.bonded[0], POLARIZE)),
		  c(parametersToBeContracted.contracted_type(i.bonded[1], POLARIZE)),
		  d(parametersToBeContracted.contracted_type(i.bonded[2], POLARIZE));
	  if ((a < ACT) && (b < ACT) && (c < ACT) && (d < ACT) /*&& (e < ACT)*/)
	  {
		  parametersToBeContracted.m_polarizes.push_back(i);
		  parametersToBeContracted.m_polarizes.back().index = a;
		  parametersToBeContracted.m_polarizes.back().bonded[0] = b;
		  parametersToBeContracted.m_polarizes.back().bonded[1] = c;
		  parametersToBeContracted.m_polarizes.back().bonded[2] = d;
	  }
  }

  for (auto const &i : m_strbends)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], STRBEND)),
      b(parametersToBeContracted.contracted_type(i.index[1], STRBEND)),
      c(parametersToBeContracted.contracted_type(i.index[2], STRBEND));
    if ((a < ACT) && (b < ACT) && (c < ACT))
    {
      parametersToBeContracted.m_strbends.push_back(i);
      parametersToBeContracted.m_strbends.back().index[0] = a;
      parametersToBeContracted.m_strbends.back().index[1] = b;
      parametersToBeContracted.m_strbends.back().index[2] = c;
    }
  }

  for (auto const &i : m_ureybrads)
  {
    std::size_t const a(parametersToBeContracted.contracted_type(i.index[0], UREYBRAD)),
      b(parametersToBeContracted.contracted_type(i.index[1], UREYBRAD)),
      c(parametersToBeContracted.contracted_type(i.index[2], UREYBRAD));
    if ((a < ACT) && (b < ACT) && (c < ACT))
    {
      parametersToBeContracted.m_ureybrads.push_back(i);
      parametersToBeContracted.m_ureybrads.back().index[0] = a;
      parametersToBeContracted.m_ureybrads.back().index[1] = b;
      parametersToBeContracted.m_ureybrads.back().index[2] = c;
    }
  }
  
  return parametersToBeContracted;
}

/*

##     ## ########  ##      ## 
##     ## ##     ## ##  ##  ## 
##     ## ##     ## ##  ##  ## 
##     ## ##     ## ##  ##  ## 
 ##   ##  ##     ## ##  ##  ## 
  ## ##   ##     ## ##  ##  ## 
   ###    ########   ###  ###  

##     ##    ###    ######## ########  ####  ######  ########  ######  
###   ###   ## ##      ##    ##     ##  ##  ##    ## ##       ##    ## 
#### ####  ##   ##     ##    ##     ##  ##  ##       ##       ##       
## ### ## ##     ##    ##    ########   ##  ##       ######    ######  
##     ## #########    ##    ##   ##    ##  ##       ##             ## 
##     ## ##     ##    ##    ##    ##   ##  ##    ## ##       ##    ## 
##     ## ##     ##    ##    ##     ## ####  ######  ########  ######  

*/

tinker::parameter::parameters::vdwc_matrices_t tinker::parameter::parameters::vdwc_matrices (void) const
{ 
  vdwc_matrices_t nbm;
  std::size_t const N(m_reduced_types.size());
  for (std::size_t i(0u); i<6u; ++i)
  {
    if (m_general.chg_scale.required(i) || m_general.vdw_scale.required(i) 
      || (m_general.vdw_scale.used(i) && i == 3U && !m_vdw14s.empty()))
    {
      double const vs = m_general.vdw_scale.factor(i), 
        cs = m_general.chg_scale.factor(i);
      nbm[i].resize(N);
      for (std::size_t row(1u), col(1u); row<N && col<N; )
      {
        combi::vdwc c(vdwc_combi(row, col, (i==3U)));
        c.E *= vs;
        c.C *= cs;
        nbm[i](row*(row+1U)/2U + col) = c;
        ++col;
        if (col>row)
        {
          col = 1u;
          ++row;
        }
      }
    }
  }
  return nbm;
}



tinker::parameter::combi::vdwc tinker::parameter::parameters::vdwc_combi (std::size_t const a, std::size_t const b, bool const use14vdw) const
{
  combi::vdwc c;
  
  if (a==b)
  {
    if (!m_charges.empty())
    {
      charge const & ca(find_chg(a));
      c.C = m_general.electric*ca.c*ca.c;
    }
    vdw const & va(use14vdw ? find_vdw14(a) : find_vdw(a));
    c.E = m_general.epsilonrule.process(va.e);
    c.R = m_general.radiusrule.process(m_general.radiussize.process(va.r));
	c.RR = va.rf;
  }
  else
  {
    if (!m_charges.empty())
    {
      charge const & ca(find_chg(a)), cb(find_chg(b));
      c.C = m_general.electric*ca.c*cb.c;
    }
    vdw const & va(use14vdw ? find_vdw14(a) : find_vdw(a)), vb(use14vdw ? find_vdw14(b) : find_vdw(b));
    c.R = m_general.radiusrule.process(m_general.radiussize.process(va.r), m_general.radiussize.process(vb.r));
    c.E = m_general.epsilonrule.process(va.e, vb.e);
	c.RR = va.rf;
  }
  if (m_general.radiustype.value == m_general.radiustype.SIGMA) c.E *= 4.0;
  return c;
}

/*

 ######  ######## ########  ########    ###    ##     ##  #######  ########   ######  
##    ##    ##    ##     ## ##         ## ##   ###   ### ##     ## ##     ## ##    ## 
##          ##    ##     ## ##        ##   ##  #### #### ##     ## ##     ## ##       
 ######     ##    ########  ######   ##     ## ## ### ## ##     ## ########   ######  
      ##    ##    ##   ##   ##       ######### ##     ## ##     ## ##              ## 
##    ##    ##    ##    ##  ##       ##     ## ##     ## ##     ## ##        ##    ## 
 ######     ##    ##     ## ######## ##     ## ##     ##  #######  ##         ######  

*/

std::ostream& tinker::parameter::operator<< (std::ostream & stream, atom const & a)
{
  stream << std::setw(12) << std::left << "atom";
  stream << std::setw(5) << std::right << a.type;
  stream << std::setw(5) << std::right << a.group;
  stream << std::setw(5) << std::right << a.symbol << "    \"";
  stream << std::setw(30) << std::left << std::string(a.description).append("\"");
  stream << std::setw(5) << std::right << a.atomic;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.mass;
  stream << std::setw(5) << std::right << a.bonds;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, ureybrad const & a)
{
  stream << std::setw(12) << std::left << "ureybrad";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  stream << std::setw(5) << std::right << a.index[2];
  stream << std::setw(20) << std::setprecision(5) << std::fixed << std::right << a.f;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.ideal;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, angle const & a)
{
  stream << std::setw(12) << std::left << "angle";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  stream << std::setw(5) << std::right << a.index[2];
  stream << std::setw(20) << std::setprecision(5) << std::fixed << std::right << a.f;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.ideal;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, bond const & a)
{
  stream << std::setw(12) << std::left << "bond";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  stream << std::setw(25) << std::setprecision(5) << std::fixed << std::right << a.f;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.ideal;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, charge const & a)
{
  stream << std::setw(12) << std::left << "charge";
  stream << std::setw(5) << std::right << a.index;
  stream << std::setw(30) << std::setprecision(5) << std::fixed << std::right << a.c;
  return stream;
}

void tinker::parameter::torsion::to_stream (std::ostream & stream, std::string type) const
{
  stream << std::setw(12) << std::left << type;
  stream << std::setw(5) << std::right << index[0];
  stream << std::setw(5) << std::right << index[1];
  stream << std::setw(5) << std::right << index[2];
  stream << std::setw(5) << std::right << index[3];
  for (std::size_t i(0u); i<number; ++i)
  {
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << force[i];
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << ideal[i];
    stream << std::setw(5) << std::right << order[i];
  }
  stream << std::endl;
}

void tinker::parameter::improper::to_stream (std::ostream & stream, std::string type) const
{
  stream << std::setw(12) << std::left << type;
  stream << std::setw(5) << std::right << center;
  stream << std::setw(5) << std::right << ligand[1];
  stream << std::setw(5) << std::right << ligand[2];
  stream << std::setw(5) << std::right << ligand[0];
  for (std::size_t i(0u); i<number; ++i)
  {
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << force[i];
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << ideal[i];
    stream << std::setw(5) << std::right << order[i];
  }
  stream << std::endl;
}

void tinker::parameter::imptor::to_stream (std::ostream & stream, std::string type) const
{
  stream << std::setw(12) << std::left << type;
  stream << std::setw(5) << std::right << ligand[0];
  stream << std::setw(5) << std::right << ligand[1];
  stream << std::setw(5) << std::right << center;
  stream << std::setw(5) << std::right << ligand[2];
  for (std::size_t i(0u); i<number; ++i)
  {
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << force[i];
    stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << ideal[i];
    stream << std::setw(5) << std::right << order[i];
  }
  stream << std::endl;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, multipole const & a)
{
  stream << std::setw(12) << std::left << "multipole";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  stream << std::setw(5) << std::right << a.index[2];
  if (a.index[3] > 0) stream << std::setw(5) << std::right << a.index[3];
  else stream << "     ";
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.charge;
  stream << std::endl;
  stream << std::setw(32) << " ";
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.dipole.x();
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.dipole.y();
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.dipole.z();
  stream << std::endl;
  stream << std::setw(32) << " ";
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(0u);
  stream << std::endl;
  stream << std::setw(32) << " ";
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(1u);
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(2u);
  stream << std::endl;
  stream << std::setw(32) << " ";
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(3u);
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(4u);
  stream << std::setw(10) << std::setprecision(5) << std::fixed << std::right << a.quadrupole(5u);
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, opbend const & a)
{
  stream << std::setw(12) << std::left << "opbend";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  std::size_t add(25);
  stream << std::setw(5) << std::right << a.index[2]; add-=5; 
  stream << std::setw(5) << std::right << a.index[3]; add-=5; 
  stream << std::setw(add) << std::setprecision(5) << std::fixed << std::right << a.f;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, polarize const & a)
{
  stream << std::setw(12) << std::left << "polarize";
  stream << std::setw(5) << std::right << a.index;
  stream << std::setw(12) << std::setprecision(5) << std::fixed << std::right << a.p;
  stream << std::setw(12) << std::setprecision(5) << std::fixed << std::right << a.pp;
  stream << std::setw(5) << std::right << a.bonded[0];
  stream << std::setw(5) << std::right << a.bonded[1];
  stream << std::setw(5) << std::right << a.bonded[2];
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, strbend const & a)
{
  stream << std::setw(12) << std::left << "strbend";
  stream << std::setw(5) << std::right << a.index[0];
  stream << std::setw(5) << std::right << a.index[1];
  stream << std::setw(5) << std::right << a.index[2];
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.f;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.ff;
  return stream;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, vdw const & a)
{
  stream << std::setw(12) << std::left << "vdw";
  stream << std::setw(5) << std::right << a.index;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.r;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.e;
  stream << std::setw(15) << std::setprecision(5) << std::fixed << std::right << a.rf;
  return stream;
}

void tinker::parameter::factors::to_stream (std::ostream & stream, std::string type) const
{
  if (order > 0) stream << type << std::left << std::setw(30-type.length()) << "-cubic" << std::fixed << std::setprecision(10) << value[0] << std::endl;
  if (order > 1) stream << type << std::left << std::setw(30-type.length()) << "-quartic" << std::fixed << std::setprecision(10) << value[1] << std::endl;
  if (order > 2) stream << type << std::left << std::setw(30-type.length()) << "-pentic" << std::fixed << std::setprecision(10) << value[2] << std::endl;
  if (order > 3) stream << type << std::left << std::setw(30-type.length()) << "-sextic" << std::fixed << std::setprecision(10) << value[3] << std::endl;
}

std::ostream& tinker::parameter::operator<< (std::ostream & stream, global const & g)
{
  stream << std::left << std::setw(30) << "forcefield" << g.forcefield << std::endl;
  stream << std::left << std::setw(30) << "vdwtype" << g.vdwtype << std::endl;
  stream << std::left << std::setw(30) << "epsilonrule" << g.epsilonrule << std::endl;
  stream << std::left << std::setw(30) << "radiusrule" << g.radiusrule << std::endl;
  stream << std::left << std::setw(30) << "radiustype" << g.radiustype << std::endl;
  stream << std::left << std::setw(30) << "radiussize" << g.radiussize << std::endl;
  g.vdw_scale.to_stream(stream, "vdw");
  g.chg_scale.to_stream(stream, "chg");
  g.mpole_scale.to_stream(stream, "mpole");
  g.polar_scale.to_stream(stream, "polar");
  g.direct_scale.to_stream(stream, "direct");
  g.mutual_scale.to_stream(stream, "mutual");
  g.bond_factor.to_stream(stream, "bond");
  g.angle_factor.to_stream(stream, "angle");
  g.opbend_factor.to_stream(stream, "opbend");
  stream << std::left << std::setw(30) << "angleindex" << g.indices[ANGLE] << std::endl;
  stream << std::left << std::setw(30) << "bondindex" << g.indices[BOND] << std::endl;
  stream << std::left << std::setw(30) << "chargeindex" << g.indices[CHARGE] << std::endl;
  stream << std::left << std::setw(30) << "improperindex" << g.indices[IMPROPER] << std::endl;
  stream << std::left << std::setw(30) << "imptorsindex" << g.indices[IMPTORS] << std::endl;
  stream << std::left << std::setw(30) << "multipoleindex" << g.indices[MULTIPOLE] << std::endl;
  stream << std::left << std::setw(30) << "opbendindex" << g.indices[OPBEND] << std::endl;
  stream << std::left << std::setw(30) << "polarizeindex" << g.indices[POLARIZE] << std::endl;
  stream << std::left << std::setw(30) << "strbendindex" << g.indices[STRBEND] << std::endl;
  stream << std::left << std::setw(30) << "torsionindex" << g.indices[TORSION] << std::endl;
  stream << std::left << std::setw(30) << "ureyindex" << g.indices[UREYBRAD] << std::endl;
  stream << std::left << std::setw(30) << "vdwindex" << g.indices[VDW] << std::endl;
  stream << std::left << std::setw(30) << "angleunit" << g.angleunit << std::endl;
  stream << std::left << std::setw(30) << "bondunit" << g.bondunit << std::endl;
  stream << std::left << std::setw(30) << "torsionunit" << g.torsionunit << std::endl;
  stream << std::left << std::setw(30) << "dielectric" << g.dielectric << std::endl;
  stream << std::left << std::setw(30) << "electric" << g.electric << std::endl;
  return stream;
}

std::ostream& tinker::parameter::combi::operator<< (std::ostream &stream, vdwc const &c)
{
  stream << std::setprecision(2) << std::scientific << "[C:" << std::setw(10) << c.C << "][E:" << std::setw(10) << c.E << "][R:" << std::setw(10) << c.R << "]";
  return stream;
}

std::ostream& tinker::parameter::combi::operator<< (std::ostream &stream, vdw const &c)
{
  stream << std::setprecision(2) << std::scientific << "[E:" << std::setw(10) << c.E << "][R:" << std::setw(10) << c.R << "]";
  return stream;
}

namespace tinker{
    namespace parameter{

        std::ostream& operator<< (std::ostream & stream, parameters const & p)
        {
          stream << p.m_general;
          for ( auto const & v : p.m_atoms ) stream << v << std::endl;
          for ( auto const & v : p.m_vdws ) stream << v << std::endl;
          for ( auto const & v : p.m_charges ) stream << v << std::endl;
          for ( auto const & v : p.m_bonds ) stream << v << std::endl;
          for ( auto const & v : p.m_angles ) stream << v << std::endl;
          for ( auto const & v : p.m_impropers ) v.to_stream(stream, "improper");
          for ( auto const & v : p.m_imptors ) v.to_stream(stream, "imptors");
          for ( auto const & v : p.m_torsions ) v.to_stream(stream, "torsion");
          for ( auto const & v : p.m_multipoles ) stream << v << std::endl;
          for ( auto const & v : p.m_opbends ) stream << v << std::endl;
          for ( auto const & v : p.m_polarizes ) stream << v << std::endl;
          for ( auto const & v : p.m_strbends ) stream << v << std::endl;
          for ( auto const & v : p.m_ureybrads ) stream << v << std::endl;
          for ( auto const & v : p.m_vdw14s ) stream << v << std::endl;
          return stream;
        }
    }

}
std::size_t tinker::parameter::parameters::req_mem (void) const
{
  //sizeof(parameters);
  std::size_t reduction_mem(m_type_of_group.capacity());
  reduction_mem += m_reduced_types.capacity();
  reduction_mem += m_reduced_groups.capacity();
  reduction_mem += contraction_map.capacity();
  reduction_mem *= sizeof(std::size_t);


  std::size_t atoms_mem(0u);
  for (auto const &i : m_atoms) atoms_mem += i.req_mem();
  std::size_t data_mem(m_vdws.capacity()*sizeof(vdw));
  data_mem += m_charges.capacity()*sizeof(charge);
  data_mem += m_bonds.capacity()*sizeof(bond);
  data_mem += m_angles.capacity()*sizeof(angle);
  data_mem += m_impropers.capacity()*sizeof(improper);
  data_mem += m_torsions.capacity()*sizeof(torsion);
  data_mem += m_imptors.capacity()*sizeof(imptor);
  data_mem += m_multipoles.capacity()*sizeof(multipole);
  data_mem += m_opbends.capacity()*sizeof(opbend);
  data_mem += m_polarizes.capacity()*sizeof(polarize);
  data_mem += m_strbends.capacity()*sizeof(strbend);
  data_mem += m_ureybrads.capacity()*sizeof(ureybrad);
  data_mem += m_vdw14s.capacity()*sizeof(m_vdw14s);
  return reduction_mem + atoms_mem + data_mem + sizeof(parameters);

}

void tinker::parameter::global::swap (global &rhs)
{
  std::swap(forcefield,rhs.forcefield);
  std::swap(vdwtype,rhs.vdwtype);
  std::swap(radiusrule,rhs.radiusrule);
  std::swap(epsilonrule,rhs.epsilonrule);
  std::swap(radiustype,rhs.radiustype);
  std::swap(radiussize,rhs.radiussize);
  std::swap(vdw_scale,rhs.vdw_scale);
  std::swap(chg_scale,rhs.chg_scale);
  std::swap(mpole_scale,rhs.mpole_scale);
  std::swap(polar_scale,rhs.polar_scale);
  std::swap(direct_scale,rhs.direct_scale);
  std::swap(mutual_scale,rhs.mutual_scale);
  std::swap(bond_factor,rhs.bond_factor);
  std::swap(angle_factor,rhs.angle_factor);
  std::swap(opbend_factor,rhs.opbend_factor);
  std::swap(angleunit,rhs.angleunit);
  std::swap(bondunit,rhs.bondunit);
  std::swap(torsionunit,rhs.torsionunit);
  std::swap(dielectric,rhs.dielectric);
  std::swap(electric,rhs.electric);
  for (std::size_t i(0u); i<KEYS; ++i) std::swap(indices[i], rhs.indices[i]);
}
