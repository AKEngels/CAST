#ifndef _GB_
#define _GB_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "configuration.h"
#include "scon_vect.h"
#include "tinker_parameters.h"
#include "coords.h"



namespace surface
{
	class surface
	{
	private:
		struct KUGEL
		{
			double   r;
			scon::vect3d<double> m;
		};

		bool   ROUND, PERMUTATE;

	public:

		std::vector<KUGEL>   K;

		// GAUSS berechnet die exakte Oberfläche von sich überschneidenden
		// Kugeln nach dem Theorem von Gauss-Bonnet.
		//
		//     literature references :
		//
		//     T.J.Richmond, "Solvent Accessible Surface Area and
		//     Excluded Volume in Proteins", Journal of Molecular Biology,
		//     178, 63 - 89 (1984)
		//
		//     L.Wesson and D.Eisenberg, "Atomic Solvation Parameters
		//     Applied to Molecular Dynamics of Proteins in Solution",
		//     Protein Science, 1, 227 - 235 (1992)
		//
		double GAUSS(scon::vect3d<double>, double, std::vector<bool>, std::vector<KUGEL>&, int, coords::Representation_3D&);
	};
}

namespace GB
{

	//
	//     ##############################################################
	//     ##                                                          ##
	//     ##  class born  --  Born radii for implicit solvation,      ##
	//     ##                  solvation energies and derivatives      ##
	//     ##                                                          ##
	//     ##############################################################
	//
	//     The methods in this class compute the Born radius of each atom,
	//     the solvation energies and derivatives for use with the various
	//     implicit solvation models.
	//
	//     literature references :
	//
	//     W.C.Still, A.Tempczyk, R.C.Hawley and T.Hendrickson,
	//     "A Semianalytical Treatment of Solvation for Molecular
	//     Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129
	//     (1990)  ("Onion" Method; see supplimentary material)
	//
	//     D.Qiu, P.S.Shenkin, F.P.Hollinger and W.C.Still, "The
	//     GB / SA Continuum Model for Solvation.A Fast Analytical Method
	//     for the Calculation of Approximate Radii", J. Phys. Chem. A,
	//     101, 3005 - 3014 (1997)  (Analytical Still Method)
	//
	//     G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized
	//     Models of Aqueous Free Energies of Solvation Based on Pairwise
	//     Descreening of Solute Atomic Charges from a Dielectric Medium",
	//     J.Phys.Chem., 100, 19824 - 19839 (1996)  (HCT Method)
	//
	//     A.Onufriev, D.Bashford and D.A.Case, "Exploring Protein
	//     Native States and Large - Scale Conformational Changes with a
	//     Modified Generalized Born Model", PROTEINS, 55, 383-394 (2004)
	//     (OBC Method)
	//
	//     T.Grycuk, "Deficiency of the Coulomb-field Approximation
	//     in the Generalized Born Model : An Improved Formula for Born
	//     Radii Evaluation", J. Chem. Phys., 119, 4817-4826 (2003)
	//     (Grycuk Method)
	//
	//     M.Schaefer, C.Bartels and M.Karplus, "Solution Conformations
	//     and Thermodynamics of Structured Peptides : Molecular Dynamics
	//     Simulation with an Implicit Solvation Model", J. Mol. Biol.,
	//     284, 835 - 848 (1998)  (ACE Method)

class born
{
private:

	struct   TINKER_prm
	{
		std::vector<double>   R_vdw;
		std::vector<size_t>   typelist;

		std::vector<tinker::parameter::bond>     bonds;
		std::vector<tinker::parameter::angle>    angles;
		std::vector<tinker::parameter::charge>   charges;
	};
	
	// Interne Parameter und Variablen ///////////////////////////
	static const coords::Coordinates  *molecule;

	static std::string   top_file, coord_file;
	static std::string   outputfile;
	static TINKER_prm    parameters;

	static size_t   N;

	static std::vector<double>   R_solv, V_solv, SA;

	static double  offset, probe, electric, off, cut;
	static double const  PI;

	// flag-Variablen
	static int   DERIV;

	enum METHODS  { VAC = -1, STILL = 0, HCT, OBC, GRYCUK, ACE, ONION, METHODNUM };
	enum SURFACES { TINKER = 0, SASASTILL, GAUSS , SURFACESNUM};
	enum RADIUS   { STD = 0, VDW };

	static METHODS    METHOD;
	static SURFACES   SURFACE;
	static RADIUS     RSOLV;

	static std::string   method, surface, rsolv, RTYPE;
	static bool          IDEAL;
	//////////////////////////////////////////////////////////////

	static scon::vect3d<double>  xyz(size_t i)
	{
		return molecule->xyz(i);
	};

	static double x(size_t i)
	{
		return xyz(i).x();
	};

	static double y(size_t i)
	{
		return xyz(i).y();
	};

	static double z(size_t i)
	{
		return xyz(i).z();
	};

	static size_t OZ(size_t i)
	{
		return molecule->atoms(i).number();
	};

	static std::vector<size_t> bonds(size_t i)
	{
		return molecule->atoms(i).bonds();
	};

	static double mass(size_t i)
	{
		return molecule->atoms(i).mass();
	};

	static double charge(size_t i)
	{
		size_t m = molecule->atoms(i).energy_type();
		
		if (m < parameters.charges.size())
		{ 
			if (parameters.charges[m].index == m) return parameters.charges[m].c;
			else if (m > 0 && (parameters.charges[m - 1].index == m)) return parameters.charges[m - 1].c;
			else if (m < parameters.charges.size() - 1 && (parameters.charges[m + 1].index == m)) return parameters.charges[m + 1].c;
		}
		for (auto &a : parameters.charges) if (a.index == m) return a.c;
	};

	// get_TINKER_prm() Liest die Gleichgewichtwerte der Bindungen,
	// Winkel und vdW-Radien aus der TINKER-prm-Datei aus
	static void get_TINKER_prm()
	{
		std::vector<tinker::parameter::angle>   angles;
		tinker::parameter::parameters           prm;

		// Löschen der alten Parameter
		parameters.R_vdw.clear();
		parameters.typelist.clear();

		prm.from_file(top_file);

		//prm.group_by_type(molecule->atoms(4).energy_type())

		// Auslesen der Typ-Liste
		for (size_t i = 0; i < N; ++i)
			parameters.typelist.push_back(prm.group_by_type(molecule->atoms(i).energy_type()));

		//for (size_t i = 0; i < N; ++i)
		//	parameters.typelist.push_back(prm.group_by_type(i+1));

		// Auslesen der Bindungslängen
		parameters.bonds = prm.bonds();

		// Auslesen der Bindungswinkel
		parameters.angles = prm.angles();

		// Auslesen der Partialladungen
		parameters.charges = prm.charges();

		// Auslesen der vdW-Radien 
		for (size_t i = 0; i < N; ++i)
			parameters.R_vdw.push_back(prm.vdws()[prm.group_by_type(molecule->atoms(i).energy_type()) - 1].r);
	}

	// sucht in den aus der TINKER-Datei eingelesenen Parametern den
	// GG-Abstand zwischen Atom a und b heraus
	static double ideal_12(size_t a, size_t b)
	{
		for (size_t k = 0; k < parameters.bonds.size(); ++k)
		{
			if (((parameters.bonds[k].index[0] == parameters.typelist[a]) && (parameters.bonds[k].index[1] == parameters.typelist[b]))
				||
				((parameters.bonds[k].index[1] == parameters.typelist[a]) && (parameters.bonds[k].index[0] == parameters.typelist[b])))
			{
				return parameters.bonds[k].ideal;
				break;
			}
		}
	}

	// sucht in den aus der TINKER-Datei eingelesenen Parametern den
	// GG-Winkel zwischen Atom a, b und c heraus
	static double ideal_13(size_t a, size_t b, size_t c)
	{
		for (size_t k = 0; k < parameters.angles.size(); ++k)
		{
			if (
				(parameters.angles[k].index[1] == parameters.typelist[b])
				&&
				(
				((parameters.angles[k].index[0] == parameters.typelist[a]) && (parameters.angles[k].index[2] == parameters.typelist[c]))
				||
				((parameters.angles[k].index[2] == parameters.typelist[a]) && (parameters.angles[k].index[0] == parameters.typelist[c]))
				)
				)
			{
				return parameters.angles[k].ideal;
				break;
			}
		}
	}

	// get_R_solv : Routine zur Bestimmung der standard Solvenzradien anhand der Ordnungszahlen
	static void get_R_solv()
	{
		size_t  oz;

		R_solv.clear();
		R_solv.resize(N);

		// Liest die Van-der-Waals Radien aus der TINKER-prm-Datei aus
		// und benutzt diese als Radien für die solvatisierten Atome
		if (RSOLV == VDW)
			R_solv = parameters.R_vdw;
		// Benutzt als Radien der solvatisierten Atome Standardradien
		else if (RSOLV == STD)
		{
			for (size_t i = 0; i < N; ++i)
			{
				oz = OZ(i);

				//Warnmeldung, falls bei einem der Atome OZ = 0 ist
				if (oz == 0) std::cout << "Warnung: Atom " << i + 1 << " wurde von CAST die Ordnungszahl 0 zugewiesen\n";

				switch (oz)
				{
				case 1:
					if (OZ(bonds(i)[0]) == 7)
						// H ist an N gebunden
						R_solv[i] = 1.15;
					else if (OZ(bonds(i)[0]) == 8)
						// H ist an O gebunden
						R_solv[i] = 1.05;
					else
						// sonst
						R_solv[i] = 1.25;
					break;
				case 3:
					R_solv[i] = 1.432;
					break;
				case 6:
					if (bonds(i).size() == 3)
						// Kohlenstoff mit 3 Bindungspartnern
						R_solv[i] = 1.875;
					else if (bonds(i).size() == 2)
						// Kohlenstoff mit 2 Bindungspartnern
						R_solv[i] = 1.825;
					else
						// sonst
						R_solv[i] = 1.9;
					break;
				case 7:
					if (bonds(i).size() == 4)
						// Stickstoff mit 4 Bindungspartnern
						R_solv[i] = 1.625;
					else if (bonds(i).size() == 1)
						// Stickstoff mit 1 Bindungspartner
						R_solv[i] = 1.6;
					else
						// sonst
						R_solv[i] = 1.7063;
					break;
				case 8:
					if (bonds(i).size() == 1)
						// Sauerstoff mit 1 Bindungspartner
						R_solv[i] = 1.48;
					else
						// sonst
						R_solv[i] = 1.535;
					break;
				case 9:
					R_solv[i] = 1.47;
					break;
				case 10:
					R_solv[i] = 1.39;
					break;
				case 11:
					R_solv[i] = 1.992;
					break;
				case 12:
					R_solv[i] = 1.7;
					break;
				case 14:
					R_solv[i] = 1.8;
					break;
				case 15:
					R_solv[i] = 1.87;
					break;
				case 16:
					R_solv[i] = 1.775;
					break;
				case 17:
					R_solv[i] = 1.735;
					break;
				case 18:
					R_solv[i] = 1.7;
					break;
				case 19:
					R_solv[i] = 2.123;
					break;
				case 20:
					R_solv[i] = 1.817;
					break;
				case 35:
					R_solv[i] = 1.9;
					break;
				case 36:
					R_solv[i] = 1.812;
					break;
				case 37:
					R_solv[i] = 2.26;
					break;
				case 53:
					R_solv[i] = 54;
					break;
				case 54:
					R_solv[i] = 1.967;
					break;
				case 55:
					R_solv[i] = 2.507;
					break;
				case 56:
					R_solv[i] = 2.188;
					break;
				default:
					R_solv[i] = 2.0;
				}
			}
		}
	}

	// "ONION"-Routinen   /////////////////////////////////////////////////////////////////
	//                                                                                   //
	// "ONION" PARAMETER und Variablen                                                   //
	static surface::surface   A;

	static void prepare_ONION()
	{
		A.K.clear();
		A.K.resize(N);

		// Radien der Atome werden eingelesen
		for (size_t i = 0; i < N; ++i)
		{
			//A.K[i].m = xyz(i);
			A.K[i].r = R_solv[i] +offset;
		}
	};

	void ONION_r()
	{
		std::vector<bool>   omit;
		double   t, fraction, area, roff, total;
		bool     done;

		alpha.clear();
		alpha.resize(N);

		omit.clear();
		omit.resize(N);
		for (size_t j = 0; j < N; ++j) omit[j] = false;

		// Mittelpunkte der Atome, die sich bei jedem Durchlauf verändern
		// werden eingelesen.
		for (size_t i = 0; i < N; ++i)
		{
			A.K[i].m = xyz(i);
		}

		for (size_t i = 0; i < A.K.size(); ++i)
		{
			t = 0.1;
			total = 0.0;
			roff = R_solv[i] + offset;

			done = false;
			omit[i] = true;
			if (i>0) omit[i - 1] = false;

			//std::cout << "Atom[" << i + 1 << "]\t";
			while (!done)
			{
				//std::cout << ".";
				roff = roff + 0.5 * t;

				area = A.GAUSS(A.K[i].m,roff,omit,A.K,0,dAES);

				fraction = area / (4.0*PI*roff*roff);
				if (fraction < 0.99)
				{
					total = total + fraction*(1.0 / (roff - 0.5*t) - 1.0 / (roff + 0.5*t));
					roff = roff + 0.5*t;
					t = 1.5*t;
				}
				else
				{
					total = total + 1.0 / (roff - 0.5*t);
					done = true;
				}
			}
			alpha[i] = 1.0 / total;
			//std::cout << "done\n";
		}
	};
	///////////////////////////////////////////////////////////////////////////////////////


	// "STILL"-Routinen ///////////////////////////////////////////////////////////////////
	//                                                                                   //
	// "STILL" Parameter und Variablen                                                   //
	std::vector<double> G_pol;

	// Atomvolumina Vi und Radien R_solv für die analytische Still-Methode werden bestimmt.
	static void prepare_STILL()
	{
		double   P_1 = 0.073,  P_2 = 0.921, P_3 = 6.211,
			     P_4 = 15.236, P_5 = 1.254;

		std::vector<size_t>    strech_12;
		double        Ri, Rj, r_ij, h_ij;

		V_solv.clear();
		V_solv.resize(N);
	

		for (size_t i = 0; i < N; ++i)
		{
			strech_12 = bonds(i);
			Ri = R_solv[i];

			V_solv[i] = (4.0 / 3.0) * PI * pow(Ri, 3.0);

			for (size_t j = 0; j < strech_12.size(); ++j)
			{
				Rj = R_solv[strech_12[j]];
				// Abstand zum gebundenen Atom wird ermittelt
				if (IDEAL)
					r_ij = 1.01 * ideal_12(i, strech_12[j]);
				else
					r_ij = 1.01 * (xyz(i) - xyz(strech_12[j])).len();

				h_ij = Ri * (1.0 + (Rj*Rj - Ri*Ri - r_ij*r_ij) / (2.0 * Ri*r_ij));

				V_solv[i] -= (PI / 3.0) * h_ij * h_ij * (3.0 * Ri - h_ij);
			}
		}
	}

	// Berechnet die Born-Radien mit den Koordinaten aus "molecule",
	// kann erst gestartet werden nachdem die Radien und Volumina
	// der Atome mindestens einmal eingelesen worden sind.
	void STILL_r()
	{
		double   P_1 = 0.073,  P_2 = 0.921, P_3 = 6.211,
			     P_4 = 15.236, P_5 = 1.254;

		std::vector<size_t>    nonbonded, strech_12, bend_13;
		double    r_ij, Ri, Rj, temp, CCF, r_a, r_b;

		G_pol.clear();
		alpha.clear();

		G_pol.resize(N);
		alpha.resize(N);
		nonbonded.resize(N);

		// Falls keine Gleichgewichtsabstände/-winkel verwendet werden,
		// müssen die Volumina, self-, 12- und 13-Terme jedes mal neu berechnet werden
		if (!IDEAL)
			prepare_STILL();

		for (size_t i = 0; i < N; ++i)
		{
			nonbonded.clear();
			nonbonded.resize(N);
			for (size_t j = 0; j < N; ++j)
				nonbonded[j] = 0;
			nonbonded[i] = 1;
			bend_13.clear();


			G_pol[i] = -166.02610865 / (R_solv[i] + offset + P_1);

			// Beiträge der an i gebundenen Atome werden berechnet
			strech_12 = bonds(i);
			for (size_t j = 0; j < strech_12.size(); ++j)
			{
				nonbonded[strech_12[j]] = 1;

				// Atome für die 1,3-Biegeterme werden ermittelt
				for (size_t k = 0; k < bonds(strech_12[j]).size(); ++k)
				{
					if (bonds(strech_12[j])[k] != i)
					{
						bend_13.push_back(strech_12[j]);
						bend_13.push_back(bonds(strech_12[j])[k]);
					}
				}

				// Abstand zum gebundenen Atom wird ermittelt
				if (IDEAL)
					r_ij = ideal_12(i, strech_12[j]);
				else
				r_ij = (xyz(i) - xyz(strech_12[j])).len();

				G_pol[i] += P_2*V_solv[strech_12[j]] / pow(r_ij, 4);
			}

			// Doppelte Einträge aus der Liste der Biegeterme werden entfernt
			// (falls 3-er oder 4-er Ringe im Molekül existieren)
			for (size_t j = 0; j < bend_13.size() / 2; ++j)
			{
				for (size_t k = j + 1; k < bend_13.size() / 2; ++k)
				{
					// Test auf 4-er Ring
					if (bend_13[j * 2 + 1] == bend_13[k * 2 + 1])
						bend_13[k * 2 + 1] = i;

					// Test auf 3-er Ring
					if ((bend_13[j * 2] == bend_13[k * 2 + 1]) && (bend_13[j * 2 + 1] == bend_13[k * 2]))
					{
						bend_13[j * 2 + 1] = i;
						bend_13[k * 2 + 1] = i;
					}
				}
			}

			// Beiträge der Biegeterme (Atome in 1,3-Stellung zu i) werden berechnet
			for (size_t j = 0; j < bend_13.size() / 2; ++j)
			{
				if (bend_13[j * 2 + 1] != i)
				{
					nonbonded[bend_13[j * 2 +1]] = 1;

					if (IDEAL)
					{
						r_a = ideal_12(i, bend_13[j * 2]);
						r_b = ideal_12(bend_13[j * 2], bend_13[j * 2 + 1]);
						r_ij = pow(r_a, 2.0) + pow(r_b, 2.0) - 2.0 *r_a*r_b*cos((2.0*PI / 360.0)*ideal_13(i, bend_13[j * 2], bend_13[j * 2 + 1]));
					}
					else
						r_ij = pow((xyz(i) - xyz(bend_13[j * 2 + 1])).len(),2.0);

					G_pol[i] += P_3*V_solv[bend_13[j * 2 + 1]] / pow(r_ij, 2.0);
				}
			}

			for (size_t j = 0; j < N; ++j)
			{
				if (nonbonded[j] == 0)
				{
					r_ij = (xyz(i) - xyz(j)).len();
					Ri = R_solv[i];
					Rj = R_solv[j];

					temp = pow(r_ij / (Ri + Rj), 2.0);
					if (temp > 1 / P_5)
						CCF = 1.0;
					else
					{
						CCF = pow(0.5*(1.0 - (cos(temp*P_5*PI))), 2.0);
					}

					G_pol[i] += P_4 * V_solv[j] * CCF / pow(r_ij, 4.0);
				}
			}
			alpha[i] = (electric/-2.0) / G_pol[i];
		}
	};
	///////////////////////////////////////////////////////////////////////////////////


	// "HCT","OBC","GRYCUK"-Routinen //////////////////////////////////////////////////
	//                                                                               //
	// "HCT"/"OBC"/"GRYCEK" Parameter und Variablen                                  //
	// S_HCT : Overlap-factor der Atome
	static std::vector<double> S_HCT, drOBC;

	// Bereitet das HCT/OBC/GRYCUK-Verfahren vor, zuvoe müssen die TINKER-Parameter eingelesen werden
	static void prepare_HCT_OBC_GRYCUK()
	{
		size_t   oz;

		// Bestimmung der overlap-Faktoren S_HCT
		S_HCT.clear();
		for (size_t i = 0; i < N; ++i)
			S_HCT.push_back(0.8);
		if (RTYPE == "MACROMODEL")
		{
			for (size_t i = 0; i < N; ++i)
			{
				oz = OZ(i);

				if (oz == 1)  S_HCT[i] = 0.85;
				if (oz == 6)  S_HCT[i] = 0.72;
				if (oz == 7)  S_HCT[i] = 0.79;
				if (oz == 8)  S_HCT[i] = 0.85;
				if (oz == 9)  S_HCT[i] = 0.88;
				if (oz == 15) S_HCT[i] = 0.86;
				if (oz == 16) S_HCT[i] = 0.96;
				if (oz == 26) S_HCT[i] = 0.88;
			}
		}

	}

	void HCT_r()
	{
		double   r, S_j, l_ij, u_ij, r_i, term , sum;

		alpha.clear();
		alpha.resize(N);


		for (size_t i = 0; i < N; ++i)
		{
			r_i = R_solv[i] + offset;
			sum = 1.0 / r_i;

			for (size_t j = 0; j < N; ++j)
			{
				if (i != j)
				{
					r = (xyz(j) - xyz(i)).len();
					S_j = S_HCT[j] * (R_solv[j] + offset);

					if (r_i < (r + S_j))
					{
						l_ij = 1.0 / std::max(r_i, std::fabs(r - S_j));
						u_ij = 1.0 / (r + S_j);
						term = l_ij - u_ij
							+ 0.25 * r * (u_ij*u_ij - l_ij*l_ij)
							+ (0.5 / r) * log(u_ij / l_ij)
							+ (0.25*S_j*S_j / r)*(l_ij*l_ij - u_ij*u_ij);

						if (r_i < (S_j - r))
						{
							term = term + 2.0 * (1.0 / r_i - l_ij);
						}

						sum = sum - 0.5 * term;
					}
				}
			}
			alpha[i] = 1.0 / sum;
		}
	};

	void OBC_r()
	{
		double   r, S_j, l_ij, u_ij, r_i, term, sum, tsum,tchain;

		alpha.clear();
		alpha.resize(N);

		if (DERIV == 1)
		{
			drOBC.clear();
			drOBC.resize(N);
			for (size_t i = 0; i < N; ++i)
			{
				drOBC[i] = 1.0;
			}
		}

		for (size_t i = 0; i < N; ++i)
		{
			r_i = R_solv[i] + offset;
			sum = 0.0;

			for (size_t j = 0; j < N; ++j)
			{
				if (i != j)
				{
					r = (xyz(j) - xyz(i)).len();
					S_j = (R_solv[j] + offset) * S_HCT[j];

					if (r_i < r + S_j)
					{
						l_ij = 1.0 / std::max(r_i, std::max(r - S_j, S_j - r));
						u_ij = 1.0 / (r + S_j);
						term = l_ij - u_ij
							+ 0.25 * r * (u_ij*u_ij - l_ij*l_ij)
							+ (0.5 / r) * log(u_ij / l_ij)
							+ (0.25*S_j*S_j / r)*(l_ij*l_ij - u_ij*u_ij);

						if (r_i < (S_j - r))
							term = term + 2.0 * (1.0 / r_i - l_ij);

						sum = sum + 0.5 * term;
					}
				}
			}
			sum = r_i * sum;
			tsum = tanh(sum - 0.8 * pow(sum, 2.0) + 4.85 * pow (sum, 3.0));
			alpha[i] = 1.0 / r_i - tsum / R_solv[i];
			alpha[i] = 1.0 / alpha[i];

			if (DERIV == 1)
			{
				tchain = r_i*(1.0 - 1.6*sum + 14.55*sum*sum);
				drOBC[i] = (1.0 - tsum*tsum)*tchain / R_solv[i];
			}
		}
	};

	// ACHTUNG: Im originalen TINKER-6.3.3.-Quelltext wird die Berechnung der Bornradien
	//          nach der Grycuk-Methode nicht durchgeführt, wenn in der keyfile SOLVATE GRYCUK
	//          steht. Der Quelltext muss abgeändert (z.B. durch Ändern der Schlüsselworte
	//          in der born-Routine in der born.f-Datei in dem man "GRYCUK" durch "OBC" und umgekehrt
	//          ersetzt) und neu compiliert werden
	void GRYCUK_r()
	{
		double   r, S_j, l_ij, u_ij, term, sum, PI34;

		alpha.clear();
		alpha.resize(N);

		PI34 = (4.0 / 3.0) * PI;

		for (size_t i = 0; i < N; ++i)
		{
			alpha[i] = 0.0;

			if (R_solv[i] > 0.0)
			{
				sum = PI34 / pow(R_solv[i], 3.0);

				for (size_t j = 0; j < N; ++j)
				{
					if ((i != j) && (R_solv[j] > 0.0))
					{
						r = (xyz(j) - xyz(i)).len();
						S_j = R_solv[j] * S_HCT[j];

						if ((R_solv[i] + r) < S_j)
						{
							l_ij = R_solv[i];
							u_ij = S_j - r;
							sum = sum + PI34 * (1.0 / pow(u_ij, 3.0) - 1.0 / pow(l_ij, 3.0));
						}
						u_ij = r + S_j;
						if ((R_solv[i] + r) < S_j)
							l_ij = S_j - r;
						else if (r < (R_solv[i] + S_j))
							l_ij = R_solv[i];
						else
							l_ij = r - S_j;

						term = (3.0*(r*r - S_j*S_j) + 6.0*u_ij*u_ij - 8.0*u_ij*r) / (r*pow(u_ij, 4.0))
							 - (3.0*(r*r - S_j*S_j) + 6.0*l_ij*l_ij - 8.0*l_ij*r) / (r*pow(l_ij, 4.0));
						sum = sum - PI*term / 12.0;
					}
				}
				// ABWEICHUNG VOM TINKER QUELLTEXT, sum kann manchmal kleiner 0 werden, wodurch für alpha NaN zurücjgegeben wird
				if (sum < 0.0) sum = sum;
				//
				
				alpha[i] = pow(sum / PI34, 1.0 / 3.0);
				if (alpha[i] <= 0.0)
					alpha[i] = 0.0001;
				alpha[i] = 1.0 / alpha[i];
			}
		}
	};

	///////////////////////////////////////////////////////////////////////////////////

	// "ACE"-Routinen   ///////////////////////////////////////////////////////////////
	//                                                                               //
	// "ACE" Parameter und Variablen /////////////
	static std::vector<double> U_ACE, W_ACE, S2_ACE;

	static void prepare_ACE()
	{
		size_t   oz, nh, m, mm;
		double   MASS;
		bool     amide;


		V_solv.clear();
		R_solv.clear();

		V_solv.resize(N);
		R_solv.resize(N);

		// Atomradien und Volumina werden bestimmt
		if (RSOLV == VDW)
			R_solv = parameters.R_vdw;
		else if (RSOLV == STD)
		{
			for (size_t i = 0; i < N; ++i)
			{
				oz = OZ(i);
				MASS = mass(i);

				if (oz == 1)
				{
					R_solv[i] = 1.468;
					V_solv[i] = 11.0;

					if ((OZ(bonds(i)[0]) == 6)
						&& (bonds(bonds(i)[0]).size() == 4))
						// Der Wasserstoff ist an einen 4-bindigen Kohlenstoff gebunden
						V_solv[i] = 11.895;
					else if ((OZ(bonds(i)[0]) == 6)
						&& (bonds(bonds(i)[0]).size() == 3))
						// Der Wasserstoff ist an einen 3-bindigen Kohlenstoff gebunden
						V_solv[i] = 13.242;
					else if ((OZ(bonds(i)[0]) == 7)
						&& (bonds(bonds(i)[0]).size() == 4))
						// Der Wasserstoff ist an einen 4-bindigen Stickstoff gebunden
					{
						V_solv[i] = 9.138;
						R_solv[i] = 0.6;
					}
					else if ((OZ(bonds(i)[0]) == 7)
						|| (OZ(bonds(i)[0]) == 8))
						// Der Wasserstoff ist an einen 3-bindigen Stickstoff oder
						// einen Sauerstoff gebunden
					{
						V_solv[i] = 9.901;
						R_solv[i] = 0.6;
					}
					else if (OZ(bonds(i)[0]) != 16)
						V_solv[i] = 13.071;
				}
				else if (oz == 6)
				{
					R_solv[i] = 2.49;
					V_solv[i] = 7.0;
					nh = 0;
					MASS = mass(i);

					// Anzahl der an Atom i gebundenen H wird ermittelt
					for (size_t j = 0; j < bonds(i).size(); ++j)
					{
						if (OZ(bonds(i)[j]) == 1)
							nh += 1;
					}

					// 4-Bindiger Kohlenstoff
					if (bonds(i).size() == 4)
					{
						if (nh == 3)
							V_solv[i] = 3.042;
						else if (nh == 2)
							V_solv[i] = 3.743;
						else if (nh == 1)
							V_solv[i] = 4.38;
					}

					// 3-Bindiger Kohlenstoff
					if (bonds(i).size() == 3)
					{
						if (nh == 1)
						{
							R_solv[i] = 2.1;
							V_solv[i] = 7.482;
						}
						else if (nh == 0)
						{
							R_solv[i] = 2.1;
							V_solv[i] = 8.288;
						}

						// ACTUNG: im TINKER-6.3.3.-Quelltext ist an dieser Stelle ein Fehler in ksolf.f
						//         - in der Subroutine "kgb" in Zeile 638 steht "k=i12(1,j)", dabei muss
						//           es "k=i12(j,i)" heißen, da k die Atomnummer des j-ten am i-ten
						//           gebundenen Atoms sein soll
						for (size_t j = 0; j < bonds(i).size(); ++j)
						{
							m = bonds(i)[j];
							if ((OZ(m) == 8) && (bonds(m).size() == 1))
							{
								R_solv[i] = 2.1;
								V_solv[i] = 7.139;
							}
						}
					}
					if ((14.5 <= MASS) && (MASS < 15.5))
					{
						R_solv[i] = 2.165;
						V_solv[i] = 33.175;
					}
					if ((13.5 <= MASS) && (MASS < 14.5))
					{
						R_solv[i] = 2.235;
						V_solv[i] = 20.862;
					}
					if ((12.5 <= MASS) && (MASS < 13.5))
					{
						R_solv[i] = 2.365;
						V_solv[i] = 11.784;
					}
				}
				else if (oz == 7)
				{
					R_solv[i] = 1.6;
					V_solv[i] = 6.0;
					nh = 0;

					// Anzahl der an Atom i gebundenen H wird ermittelt
					for (size_t j = 0; j < bonds(i).size(); ++j)
					{
						if (OZ(bonds(i)[j]) == 1)
							nh = nh + 1;
					}

					if (bonds(i).size() == 4)
					{
						if (nh == 3)
							V_solv[i] = 2.549;
						else if (nh == 2)
							V_solv[i] = 3.304;
					}
					else if (bonds(i).size() == 3)
					{
						amide = false;
						for (size_t j = 0; j < bonds(i).size(); ++j)
						{
							m = bonds(i)[j];
							if (OZ(m) == 6)
							{
								for (size_t k = 0; k < bonds(m).size(); ++k)
								{
									mm = bonds(m)[k];
									if ((OZ(mm) == 8)
										&& (bonds(mm).size() == 1)
										&& (i != mm))
										amide = true;
								}
							}
						}
						if (amide)
						{
							if (nh == 0) V_solv[i] = 7.189;
							else if (nh == 1) V_solv[i] = 6.03;
							else if (nh == 2) V_solv[i] = 5.693;
						}
						else
						{
							// An dieser Stelle ist im TINKER-6.3.3.-Quelltext ein Fehler
							// sowohl im if als auch im else Argument stehr (nh == 2)
							if (nh == 2) V_solv[i] = 5.677;
							else if (nh != 2) V_solv[i] = 6.498;
						}
					}
				}
				else if (oz == 8)
				{
					R_solv[i] = 1.6;
					V_solv[i] = 12.0;

					if (bonds(i).size() == 1)
					{
						V_solv[i] = 13.532;
						m = bonds(i)[0];

						if (OZ(m) == 15)
							V_solv[i] = 17.202;
						else
						{
							for (size_t j = 0; j < bonds(m).size(); ++j)
							{
								mm = bonds(m)[j];

								if ((OZ(mm) == 8)
									&& (bonds(mm).size() == 1)
									&& (i != mm))
									V_solv[i] = 15.4;
							}
						}
					}
					if (bonds(i).size() == 2)
					{
						V_solv[i] = 10.642;
						for (size_t j = 0; j < bonds(i).size(); ++j)
						{
							m = bonds(i)[j];
							if (OZ(m) == 15)
								V_solv[m] = 11.416;
						}
					}
				}
				else if (oz == 12)
				{
					R_solv[i] = 1.0;
					V_solv[i] = 15.235;
				}
				else if (oz == 15)
				{
					R_solv[i] = 1.89;
					V_solv[i] = 6.131;
				}
				else if (oz == 16)
				{
					R_solv[i] = 1.89;
					V_solv[i] = 17.232;
					for (size_t j = 0; j < bonds(i).size(); ++j)
					{
						m = bonds(i)[j];
						if (OZ(m) == 16) V_solv[i] = 18.465;
					}
				}
				else if (oz == 26)
				{
					R_solv[i] = 0.65;
					V_solv[i] = 9.951;
				}
				else
				{
					// Abweichend vom TINKER-Quelltext werden, falls keine passenden Atomtypen gefungen wurden,
					// die Radien und Volumina nicht auf 0.0 gesetzt, da dies zu einer Division durch Null führt,
					R_solv[i] = 2.0;
					V_solv[i] = 8.0;
				}
			}
		}

		// Bestimmung der paarweisen Parameter der ACE-Methode
		double   c1, c2, c3, pi2, ri, rj, vj, a, tij, fij, qij, qterm, omgij, s2ij, s3ij;
		size_t   ic, jc, nc;

		c1 = 4.0 / (3.0 * PI);
		c2 = 77.0 * PI * pow(2.0, 0.5) / 512.0;
		c3 = 2.0 * PI * pow(PI, 0.5);
		pi2 = 1.0 / (PI * PI);

		W_ACE.clear();
		U_ACE.clear();
		S2_ACE.clear();

		nc = parameters.typelist.size();
		W_ACE.resize(nc*nc);
		U_ACE.resize(nc*nc);
		S2_ACE.resize(nc*nc);

		for (size_t i = 0; i < N; ++i)
		{
			ri = R_solv[i];
			ic = parameters.typelist[i];

			for (size_t j = 0; j < N; ++j)
			{
				jc = parameters.typelist[j];
				// ACTUNG: Im TINKER-6.3.3.-Quelltext ist hier ein Fehler in der ksolv.f Routine
				//         "kgb" (Zeile 758 und 759) , es steht die Klasse des Atoms und nicht seine Nummer
				//         im Argument von "rsolv(kc)" und "vsolv(kc)", dabei muss es "rsolv(k)" und "vsolv(k)" heißen
				rj = R_solv[j];
				vj = V_solv[j];
				a = std::max(1.2, ri / rj);

				tij = 0.5 * PI *(a*a * rj*rj / (ri*ri));
				fij = 2.0 / (1.0 + tij) - 1.0 / (1.0 + 2.0*tij);
				qij = tij * pow(1.0 / (1.0 + 2.0*tij), 0.5);
				qterm = qij - atan(qij);

				if (j != i)
					omgij = vj * qterm * pi2 / pow(a*rj, 4.0);
				else
					omgij = c1 * qterm / (ri * pow(a, 4.0));

				s2ij = 3.0 * qterm * pow(a*rj, 2.0) / ((3.0 + fij) * qij - 4.0*atan(qij));
				s3ij = s2ij * pow(s2ij, 0.5);
				
				U_ACE[(ic - 1)*nc + jc -1] = c2 * ri / (1.0 - (c3*s3ij*ri*omgij / vj));
				W_ACE[(ic - 1)*nc + jc -1] = omgij;
				S2_ACE[(ic - 1)*nc + jc -1] = s2ij;
			}
		}

	};

	void ACE_r()
	{
		double   r, r_i, term, gself, b0, expterm, rmu;
		size_t   it, jt, nc;

		nc = parameters.typelist.size();
		alpha.clear();
		alpha.resize(N);

		b0 = 0.0;
		for (size_t i = 0; i < N; ++i)
			b0 += V_solv[i];

		b0 = pow(0.75*b0 / PI, 1.0 / 3.0);
		for (size_t i = 0; i < N; ++i)
		{
			r_i = R_solv[i];
			it = parameters.typelist[i];

			gself = 1.0 / r_i + 2.0 * W_ACE[(it - 1)*nc + it - 1];
			for (size_t j = 0; j < N; ++j)
			{
				if (j != i)
				{
					r = (xyz(j) - xyz(i)).len();
					jt = parameters.typelist[j];

					expterm = W_ACE[(it - 1)*nc + jt - 1] * exp(-1.0 * r*r / S2_ACE[(it - 1)*nc + jt - 1]);
					rmu = pow(r, 4.0) + pow(U_ACE[(it - 1)*nc + jt - 1], 4.0);
					term = (V_solv[j] / (8.0*PI)) * pow(r*r*r / rmu, 4.0);
					gself = gself - 2.0 * (expterm + term);
				}
			}
			if (gself >= 0.5 / b0)
				alpha[i] = 1.0 / gself;
			else
				alpha[i] = 2.0 * b0 * (1.0 + b0*gself);
		}
	};


	//////////////////////////////////////////////////////////////////////////////////

	// Routinen zur Berechnung der freien Lösungsenthalpie ///////////////////////////
	static double   A_solv;

	static void prepare_E()
	{
		// Oberflächenfaktoren für den nicht-polaren Energieterm
		switch (METHOD)
		{
		case ONION:
			A_solv = 0.0072;
			break;
		case STILL:
			A_solv = 0.0049;
			break;
		case HCT:
		case OBC:
		case GRYCUK:
			A_solv = 0.0054;
			break;
		case ACE:
			A_solv = 0.003;
		}
	};

	// Berechnung der polaren Energieterme
	void get_G_pol()
	{
		double c0 =  3.76033795442375E2  ;
		double c1 = -2.41328018087701E-9 ;
		double c2 =  6.12601892068780E-21;
		double c3 = -7.65971363122508E-33;
		double c4 =  4.71232224668292E-45;
		double c5 = -1.14238115071101E-57;

		double f0 = -9.55179508356813E-10;
		double f1 =  7.89189860088975E-21;
		double f2 = -2.75785243092308E-32;
		double f3 =  5.28011237195066E-44;
		double f4 = -5.97732747832433E-56;
		double f5 =  3.99786257871453E-68;
		double f6 = -1.46144156032738E-80;
		double f7 =  2.24975609656308E-93;

		double   dwater = 78.3, f = -electric*(1.0-1.0/dwater), 
			     r, r2, r3, r4, r5, r6, r7, rb, fik, fgb, fgb2,
				 e, rm2, fgm, shift, taper, trans, expterm, de,
				 derb, dtaper, dtrans, xr, yr, zr, dedx, dedy, dedz;

		for (size_t i = 0; i < N; ++i)
		{
			if (charge(i) != 0.0)
			{
				for (size_t k = i; k < N; ++k)
				if (charge(k) != 0.0)
				{
					interactions++;

					r = (xyz(i) - xyz(k)).len();
					r2 = r*r;

					if (r <= off)
					{
						fik = f * charge(i) * charge(k);
						rb = alpha[i] * alpha[k];
						expterm = exp(-0.25*r2 / rb);
						fgb2 = r2 + rb*expterm;
						fgb = pow(fgb2, 0.5);
						e = fik / fgb;

						if (DERIV == 1)
						{
							de = -e*(r - 0.25*r*expterm) / fgb2;
							derb = -e*expterm*(0.5 + 0.125*r2 / rb) / fgb2;
						}

						rm2 = pow(0.5*(off + cut), 2.0);
						fgm = pow(rm2 + rb*exp(-0.25*rm2 / rb), 0.5);
						shift = fik / fgm;
						e = e - shift;
						if (r > cut)
						{
							r3 = r*r2;
							r4 = r2*r2;
							r5 = r2*r3;
							r6 = r3*r3;
							r7 = r3*r4;

							taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0;
							trans = fik*(f7*r7 + f6*r6 + f5*r5 + f4*r4 + f3*r3 + f2*r2 + f1*r + f0);
							e = e*taper + trans;

							if (DERIV == 1)
							{
								dtaper = 5.0*c5*r4 + 4.0*c4*r3 + 3.0*c3*r2 + 2.0*c2*r + c1;
								dtrans = fik*(7.0*f7*r6 + 6.0*f6*r5 + 5.0*f5*r4 + 4.0*f4*r3 + 3.0*f3*r2 + 2.0*f2*r + f1);
								derb = derb*taper;
								de = e*dtaper + de*taper + dtrans;
							}
						}
						if (i == k)
						{
							ES = ES + 0.5*e;
							AES[i] = AES[i] + 0.5*e;

							if (DERIV == 1)
							{
								dalpha[i] = dalpha[i] + derb * alpha[i];
							}
						}
						else
						{
							ES = ES + e;
							AES[i] = AES[i] + 0.5*e;
							AES[k] = AES[k] + 0.5*e;

							if (DERIV == 1)
							{
								de = de / r;

								xr = x(i) - x(k);
								yr = y(i) - y(k);
								zr = z(i) - z(k);

								dedx = de * xr;
								dedy = de * yr;
								dedz = de * zr;

								update_dAES(i,k,dedx,dedy,dedz);

								dalpha[i] = dalpha[i] + derb * alpha[k];
								dalpha[k] = dalpha[k] + derb * alpha[i];

								update_vir(dedx,dedy,dedz,xr,yr,zr);
							}
						}
					}
				}
			}
		}
	}

	void get_AES()
	{
		double   term = 4.0*PI, e;

		AES.clear();
		AES.resize(N);
		ES = 0.0;
		interactions = N;

		if (DERIV == 1)
		{
			dalpha.clear();
			dalpha.resize(N);
			dAES.assign(N,coords::Cartesian_Point());

			for (auto &x : vir) for (auto &y : x) y = 0.0;
		}

		// Born-Radien werden berechnet
		get_r();
		
		// Lösungsmittelzugängliche Oberfläche wird berechnet
		switch (SURFACE)
		{
		case TINKER:
			// Bornradienabhängig
			get_SA_TINKER();
			break;
		case SASASTILL:
			get_SA_STILL();
			break;
		case GAUSS:
			get_SA_GAUSS();
		}

		// Der unpolare Energieterm wird berechnet
		for (size_t i = 0; i < N; ++i)
		{
			e = A_solv * SA[i];

			AES[i] += e;
			ES += e;

			if (DERIV == 1)
			{
				if ((SURFACE == SASASTILL) || (SURFACE == GAUSS)) dAES[i] *= A_solv;
				if (alpha[i] != 0.0) dalpha[i] = dalpha[i] - 6.0*e / alpha[i];
			}
		}

		// Berechnung des polaren Energieterms
		get_G_pol();

		// Ketten-Regel Komponeneten (Ableitungen der Bornradien) der polaren Energieterme werden berechnet
		if (DERIV == 1)
		{
			switch (METHOD)
			{
			case STILL:
				STILL_dr();
				break;
			case HCT:
			case OBC:
				//if ((METHOD == OBC) && (SURFACE != TINKER)) prepare_HCT_OBC_GRYCUK();
				HCT_OBC_dr();
				break;
			case GRYCUK:
				GRYCUK_dr();
				break;
			case ONION:
				if (SURFACE != GAUSS)
				{
					prepare_STILL();
					STILL_dr();
				}
			}
		}
	};

	// Born-Radien Ketten-Regel Komponeneten  /////////////////////////////////////////////////

	// Ketten-Regel-Komponenten der STILL-Methode                                           ///
	void STILL_dr()
	{
		double   P_1 = 0.073,  P_2 = 0.921, P_3 = 6.211,
			     P_4 = 15.236, P_5 = 1.254;

		double   gpi, r, theta, cosq, term, ccf, dccf,
			     ratio, de, xr, yr, zr, dedx, dedy, dedz;
		std::vector<double>   skip(N);

		for (size_t i = 0; i < N; ++i)
			skip[i] = -1;

		for (size_t i = 0; i < N; ++i)
		{
			// Atome für nichtbindende Wechselwirkungen erden ermittelt
			for (size_t j = 0; j < bonds(i).size(); ++j)
			{
				skip[bonds(i)[j]] = i;
				for (size_t k = 0; k < bonds(bonds(i)[j]).size(); ++k)
					skip[bonds(bonds(i)[j])[k]] = i;
			}
			gpi = 2.0*alpha[i] * alpha[i] / electric;

			for (size_t k = 0; k < N; ++k)
			{
				if (skip[k] != i)
				{
					r = (xyz(k) - xyz(i)).len();

					ratio = r*r / pow(R_solv[i] + R_solv[k], 2.0);
					if (ratio > 1.0 / P_5)
					{
						ccf = 1.0;
						dccf = 0.0;
					}
					else
					{
						theta = ratio * PI * P_5;
						cosq = cos(theta);
						term = 0.5*(1.0 - cosq);
						ccf = term*term;
						dccf = 2.0*term*sin(theta)*PI*P_5*ratio;
					}
					de = dalpha[i] * P_4*gpi*V_solv[k] * (4.0*ccf - dccf) / pow(r, 6.0);

					xr = x(k) - x(i);
					yr = y(k) - y(i);
					zr = z(k) - z(i);

					dedx = de * xr;
					dedy = de * yr;
					dedz = de * zr;

					update_dAES(i, k, dedx, dedy, dedz);
					update_vir(dedx, dedy, dedz, xr, yr, zr);
				}
			}
		}

	};

	// Ketten-Regel-Komponenten der HCT/OBC-Methode                                           ///
	void HCT_OBC_dr()
	{
		double   lik, lik2, lik3, uik, uik2, uik3, dlik, duik,
			     r, ri, rk, sk, t1, t2, t3, de, xr, yr, zr, rb2,
				 dedx, dedy, dedz;

		for (size_t i = 0; i < N; ++i)
		{
			ri = R_solv[i] + offset;
			for (size_t k = 0; k < N; ++k)
			{
				if (i != k)
				{
					rk = R_solv[k] + offset;
					sk = S_HCT[k] * rk;
					r = (xyz(k) - xyz(i)).len();

					if (ri < r + sk)
					{
						lik = 1.0 / std::max(ri, std::fabs(r - sk));
						uik = 1.0 / (r + sk);

						lik2 = lik*lik;
						uik2 = uik*uik;
						lik3 = lik*lik2;
						uik3 = uik*uik2;

						dlik = 1.0;
						if (ri >= r - sk) dlik = 0.0;

						duik = 1.0;

						t1 = 0.5*lik2 + 0.25*sk*sk*lik3 / r - 0.25*(lik / r + lik3*r);
						t2 = -0.5*uik2 - 0.25*sk*sk*uik3 / r + 0.25*(uik / r + uik3*r);
						t3 = 0.125*(1.0 + sk*sk / (r*r))*(lik2 - uik2) + 0.25*log(uik / lik) / (r*r);

						rb2 = alpha[i] * alpha[i];

						if (METHOD == OBC) rb2 *= drOBC[i];

						de = dalpha[i] * rb2 * (dlik*t1 + duik*t2 + t3) / r;

						xr = x(k) - x(i);
						yr = y(k) - y(i);
						zr = z(k) - z(i);

						dedx = de * xr;
						dedy = de * yr;
						dedz = de * zr;

						update_dAES(i, k, dedx, dedy, dedz);
						update_vir(dedx, dedy, dedz, xr, yr, zr);
					}
				}
			}
		}
	};

	void GRYCUK_dr()
	{
		double   pi43 = 4.0*PI / 3.0, factor = -pow(PI, 1.0 / 3.0)*pow(6.0, 2.0 / 3.0) / 9.0,
			     term, r, sk, de, uik, lik, drb, dedx, dedy, dedz, xr, yr, zr;

		for (size_t i = 0; i<N; ++i)
		{
			if (R_solv[i] > 0.0)
			{
				term = pi43 / pow(alpha[i], 3.0);
				term = factor / pow(term, 4.0 / 3.0);

				for (size_t k = 0; k < N; ++k)
				{
					if ((i != k) && (R_solv[k] > 0.0))
					{
						sk = R_solv[k] * S_HCT[k];
						r = (xyz(k) - xyz(i)).len();
						de = 0.0;

						if (R_solv[i] + r < sk)
						{
							uik = sk - r;
							de = -4.0*PI / pow(uik, 4.0);
						}
						// ACTUNG: Fehler im TINKER-4.3.0-Quelltext, das Argument der
						//         zweiten if-Anweisung entspricht dem der vorhergehenden
						if (R_solv[i] + r < sk)
						{
							lik = sk - r;
							de = de + 0.25*PI*(sk*sk - 4.0*sk*r + 17.0*r*r) / (r*r*pow(lik, 4.0));
						}
						else if (r < R_solv[i] + sk)
						{
							lik = R_solv[i];
							de = de + 0.25*PI*(2.0*R_solv[i] * R_solv[i] - sk*sk - r*r) / (r*r*pow(lik, 4.0));
						}
						else
						{
							lik = r - sk;
							de = de + 0.25*PI*(sk*sk - 4.0*sk*r + r*r) / (r*r*pow(lik, 4.0));
						}
						uik = r + sk;
						de = de - 0.25*PI*(sk*sk + 4.0*sk*r + r*r) / (r*r*pow(uik, 4.0));
						drb = term*de / r;
						de = drb*dalpha[i];

						xr = x(k) - x(i);
						yr = y(k) - y(i);
						zr = z(k) - z(i);

						dedx = de*xr;
						dedy = de*yr;
						dedz = de*zr;

						update_dAES(i, k, dedx, dedy, dedz);
						update_vir(dedx, dedy, dedz, xr, yr, zr);
					}
				}
			}
		}
	};

	void ACE_dr()
	{
		size_t   it, kt, nc = parameters.typelist.size();
		double   rbi2, r, r2, r3, r4, r6, xr, yr, zr,
			     dedx, dedy, dedz, s2ik, ws2, uik4,
				 ratio, rusum, expterm, de, de1, de2;

		for (size_t i = 0; i < N; ++i)
		{
			it = parameters.typelist[i];
			rbi2 = alpha[i] * alpha[i];
			for (size_t k = 0; k < N; ++k)
			{
				if (i != k)
				{
					kt = parameters.typelist[k];
					s2ik = 1.0 / S2_ACE[(it - 1)*nc + kt - 1];
					ws2 = W_ACE[(it - 1)*nc + kt - 1] * s2ik;
					uik4 = pow(U_ACE[(it - 1)*nc + kt - 1], 4.0);
					r = (xyz(i) - xyz(k)).len();
					r2 = r*r;
					r3 = r*r2;
					r4 = r2*r2;
					r6 = r3*r3;

					rusum = r4 + uik4;
					ratio = r3 / rusum;
					expterm = exp(-r2*s2ik);
					de1 = -4.0*r*ws2*expterm;
					de2 = 3.0*r2 / rusum - 4.0*r6 / pow(rusum, 2.0);
					de = dalpha[i] * alpha[i] * alpha[i] * (de1 + V_solv[k] * pow(ratio, 3.0)*de2 / PI) / r;

					xr = x(i) - x(k);
					yr = y(i) - y(k);
					zr = z(i) - z(k);

					dedx = de*xr;
					dedy = de*yr;
					dedz = de*zr;

					update_dAES(i, k, dedx, dedy, dedz);
					update_vir(dedx, dedy, dedz, xr, yr, zr);
				}
			}
		}
	};

	void update_dAES(size_t &i, size_t &k, double &dedx, double &dedy, double &dedz)
	{
		dAES[i].setX(dAES[i].x() + dedx);
		dAES[i].setY(dAES[i].y() + dedy);
		dAES[i].setZ(dAES[i].z() + dedz);

		dAES[k].setX(dAES[k].x() - dedx);
		dAES[k].setY(dAES[k].y() - dedy);
		dAES[k].setZ(dAES[k].z() - dedz);
	};

	void update_vir(double &dedx, double &dedy, double &dedz, double &xr, double &yr, double zr)
	{
		vir[0][0] = vir[0][0] + xr * dedx;
		vir[1][0] = vir[1][0] + yr * dedx;
		vir[2][0] = vir[2][0] + zr * dedx;
		vir[0][1] = vir[0][1] + xr * dedy;
		vir[1][1] = vir[1][1] + yr * dedy;
		vir[2][1] = vir[2][1] + zr * dedy;
		vir[0][2] = vir[0][2] + xr * dedz;
		vir[1][2] = vir[1][2] + yr * dedz;
		vir[2][2] = vir[2][2] + zr * dedz;
	};
	
	// Methoden zur berechnung der lösungsmittelzugänglichen Oberfläche der Atome
	void get_SA_TINKER(void);
	void get_SA_STILL (void);
	void get_SA_GAUSS (void);

	/////////////////////////////////////////////////////////////////////////////////

	static void clean_methods();

	static void prepare()
	{
		// Parameter aus der TINKER-prm DAtei werden eingelesen
		get_TINKER_prm();

		clean_methods();
		// Methoden zur Berechnung der Born-Radien werden vorbereitet
		switch (METHOD)
		{
		case STILL:
			get_R_solv();
			prepare_STILL();
			break;
		case HCT:
		case OBC:
		case GRYCUK:
			get_R_solv();
			prepare_HCT_OBC_GRYCUK();
			break;
		case ACE:
			prepare_ACE();
			break;
		case ONION:
			get_R_solv();
			prepare_ONION();
		}

		prepare_E();
	};

public:

	// Variablenblock   //////////////////////////////////
	//
	// Vektor für die Born-Radien
	std::vector<double>   alpha;
	// Ableitungen der Bornradien
	std::vector<double>   dalpha;
	// Vektor für die Lösungsenergien
	std::vector<double>   AES;
	// Ableitungen der Lösungsenergien
	//std::vector<double>   dAES_x, dAES_y, dAES_z;
	coords::Representation_3D   dAES;
	// Virial Tensor
	std::array<std::array<double, 3>, 3>   vir;

	double ES;
	long   interactions;
	//////////////////////////////////////////////////////

	born(int DERIV)
	{
		set_DERIV(DERIV);
	};

	born(coords::Coordinates &coords)
	{
		molecule = &coords;
		N = molecule->size();
		
		top_file   = Config::get().general.paramFilename;
		coord_file = Config::get().general.inputFilename;

		DERIV = 0;

		prepare();
	};

	static void set(coords::Coordinates &coords)
	{
		molecule = &coords;
		N = molecule->size();

		top_file = Config::get().general.paramFilename;
		coord_file = Config::get().general.inputFilename;

		DERIV = 0;

		prepare();
	};

	static void set_METHOD(METHODS new_METHOD)
	{
		switch (new_METHOD)
		{
		case STILL:
			METHOD = STILL;
			method = "STILL analytisch";
			break;
		case HCT:
			METHOD = HCT;
			method = "HCT";
			break;
		case OBC:
			METHOD = OBC;
			method = "OBC";
			break;
		case GRYCUK:
			METHOD = GRYCUK;
			method = "GRYCUK";
			break;
		case ACE:
			METHOD = ACE;
			method = "ACE";
			break;
		case ONION:
			METHOD = ONION;
			method = "STILL numerisch";
			break;
		default:
			METHOD = STILL;
		}
		prepare();
	};

	static void SET_METHOD()
	{
		switch (Config::get().general.solvationmethod)
		{
		case config::solvs::STILL:
			set_METHOD(METHODS::STILL);
			break;
		case config::solvs::HCT:
			set_METHOD(METHODS::HCT);
			break;
		case config::solvs::OBC:
			set_METHOD(OBC);
			break;
		case config::solvs::GRYCUK:
			set_METHOD(METHODS::GRYCUK);
			break;
		case config::solvs::ACE:
			set_METHOD(METHODS::ACE);
			break;
		case config::solvs::ONION:
			set_METHOD(METHODS::ONION);
			break;
		default:
			set_METHOD(METHODS::STILL);
		}
	};

	static bool set_method(std::string new_method)
	{
		if (new_method == "ONION")
		{
			set_METHOD(ONION);
			return true;
		}
		else if (new_method == "STILL")
		{
			set_METHOD(STILL);
			return true;
		}
		else if (new_method == "HCT")
		{
			set_METHOD(HCT);
			return true;
		}
		else if (new_method == "OBC")
		{
			set_METHOD(OBC);
			return true;
		}
		else if (new_method == "GRYCUK")
		{
			set_METHOD(GRYCUK);
			return true;
		}
		else if (new_method == "ACE")
		{
			set_METHOD(ACE);
			return true;
		}
		else
			return false;
	};



	static void set_ideal(bool value)
	{
		IDEAL = value;
		prepare();
	};

	static void set_DERIV(int D)
	{
		if (D == 0) DERIV = 0;
		if (D == 1) DERIV = 1;
	};

	static bool set_rvdw(std::string value)
	{
		if (value == "STD")
		{
			RSOLV = STD;
			rsolv = "TINKER standard vdW-Radien";
			prepare();
			return true;
		}
		else if (value == "VDW")
		{
			RSOLV = VDW;
			rsolv = "vdW-Radien der Parameterdatei";
			prepare();
			return true;
		}
		else return false;
	};

	static bool set_surface(std::string value)
	{
		if (value == "TINKER")
		{
			SURFACE = TINKER;
			surface = "TINKER-Methode";
			return true;
		}
		else if (value == "STILL")
		{
			SURFACE = SASASTILL;
			surface = "Still-Methode";
			return true;
		}
		else if (value == "GAUSS")
		{
			SURFACE = GAUSS;
			surface = "Theorem von Gauss-Bonnet";
			return true;
		}
		else return false;
	};

	static void SET_SURFACE()
	{
		switch (Config::get().general.surfacemethod)
		{
		case config::surfs::TINKER:
			GB::born::set_surface("TINKER");
			break;
		case config::surfs::SASASTILL:
			GB::born::set_surface("STILL");
			break;
		case config::surfs::GAUSS:
			GB::born::set_surface("GAUSS");
			break;
		}
	};

	// Berechnet die Born-Radien nach der ausgewählten Methode
	std::vector<double> get_r()
	{
		switch (METHOD)
		{
		case STILL:
			STILL_r();
			break;
		case HCT:
			HCT_r();
			break;
		case OBC:
			OBC_r();
			break;
		case GRYCUK:
			GRYCUK_r();
			break;
		case ACE:
			ACE_r();
			break;
		case ONION:
			ONION_r();
		}

		for (size_t i = 0; i < alpha.size(); ++i)
		{
			if ((alpha[i]<0.0) || (alpha[i]>500.0)) alpha[i] = 500.0;
		}

		write_file("test_r.rtf");

		return alpha;
	};

	double get_E()
	{
		get_AES();

		write_file("test.rtf");

		return ES;
	};
	

	// Ein-/Ausgabe Routinen   /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void write()
	{
		std::cout << '\n' << "   Born-Radien der Atome nach der \""<< method<<"\"-Methode [A]:" << '\n' << '\n';
		for (size_t i = 0; i < alpha.size(); ++i)
			std::cout << "      a_" << std::setw(6) << std::left << i + 1 << " = " << std::setw(20) << std::right << alpha[i] << '\n';
	};

	void write_ES()
	{
		std::cout << '\n' << "   Solvatisierungsenergie:   ES = " << ES<<'\n'<<'\n';
		for (size_t i = 0; i < AES.size(); ++i)
		{
			std::cout << "      AES[" << i + 1 << "] \t= " << AES[i] << '\n';
		}

		if (DERIV == 1)
		{
			std::cout <<'\n'<< "   Ableitungen der Bornradien:" << '\n' << '\n';
			for (size_t i = 0; i < N; ++i)
				std::cout << "      drbb(" << i + 1 << ")\t = " << dalpha[i] << '\n';
			std::cout <<'\n'<< "   Ableitungen der Energien der Atome nach kartesischen Koordinaten:" << '\n' << '\n';
			for (size_t i = 0; i < N; ++i)
			{
				std::cout << "      dAES_x(" << i + 1 << ")= " << dAES[i].x() <<
					"   \t dAES_y(" << i + 1 << ")= " << dAES[i].y() <<
					"   \t dAES_z(" << i + 1 << ")= " << dAES[i].z() << '\n';
			}
		}
	};

	void write_file(std::string outputfilename)
	{
		std::ofstream   out;

		out.open(outputfilename);

		out << '\n'
			<< "GB/SA" << '\n' << '\n'
			<< "   Methode               : " << method << '\n'
			<< "   Oberflächenberechnung : " << surface << '\n'
			<< "   Radien der Atome      : " << rsolv << '\n' << '\n'
			<< "   Topologiedatei   : " << top_file << '\n'
			<< "   Koordinatendatei : " << coord_file << '\n' << '\n';

		out << "Radien der solvatisierten Atome [A]:" << '\n' << '\n';
		for (size_t i = 0; i < R_solv.size(); ++i)
		{
			out << "   " << "r_" << std::setw(7) << std::left << i + 1 << "= " << std::setw(20) << std::right << std::fixed << std::setprecision(20) << R_solv[i] << '\n';
		}

		if ((METHOD == STILL) || (METHOD == ACE))
		{
			out << '\n' << "Volumina der solvatisierten Atome [A^3]:" << '\n' << '\n';
			for (size_t i = 0; i < V_solv.size(); ++i)
			{
				out << "   " << "V_" << std::setw(7) << std::left << i + 1 << "= " << std::setw(20) << std::right << std::fixed << std::setprecision(20) << V_solv[i] << '\n';
			}
		}

		if ((method == "HCT") || (method == "OBC") || (method == "GRYCUK"))
		{
			for (size_t i = 0; i < S_HCT.size(); ++i)
				out << "   " << "SHCT(" << std::setw(7) << std::left << i + 1 << ") \t= " << std::setw(20) << std::right << std::fixed << std::setprecision(20) << S_HCT[i] << '\n';
		}

		out << '\n' << "Born-Radien der solvatisierten Atome [A]:" << '\n' << '\n';
		for (size_t i = 0; i < alpha.size(); ++i)
		{
			out << "   " << "a_" << std::setw(7) << std::left << i + 1 << "= " << std::setw(20) << std::right << std::fixed << std::setprecision(20) << alpha[i] << '\n';
		}

		out << '\n' << "Solvatisierungsenergie:   ES = " << ES << '\n'<<'\n';
		out << '\n' << "Solvatisierungsenergie der einzelnen Atome:" << '\n' << '\n';
		for (size_t i = 0; i < AES.size(); ++i)
		{
			out << "   AES[" << i + 1 << "]  \t= " << AES[i] << '\n';
		}

		if (DERIV == 1)
		{
			out << '\n' << "Ableitungen der Bornradien:" << '\n' << '\n';
			for (size_t i = 0; i < N; ++i)
				out << "   drbb(" << i + 1 << ")   \t = " << dalpha[i] << '\n';
			out << '\n' << "Ableitungen der Energien der Atome nach kartesischen Koordinaten:" << '\n' << '\n';
			for (size_t i = 0; i < N; ++i)
			{
				out << "   dAES_x(" << i + 1 << ")  \t= " << dAES[i].x() <<
					"   \t dAES_y(" << i + 1 << ")  \t= " << dAES[i].y() <<
					"   \t dAES_z(" << i + 1 << ")  \t= " << dAES[i].z() << '\n';
			}
		}

		out.close();
	};

	

};
}

#endif