#include "configuration.h"
//#include "gbsa.h"
//
//namespace GB
//{
//	// Initialisierung von static-Variablen   //////////////////////////////
//	const coords::Coordinates* born::molecule = NULL;
//	size_t born::N = 0;
//
//	born::TINKER_prm born::parameters;
//
//	std::string   born::top_file   = "";
//	std::string   born::coord_file = "";
//	std::string   born::outputfile = "abl.rtf";
//
//	std::vector<double> born::R_solv;
//	std::vector<double> born::V_solv;
//	std::vector<double> born::SA;
//	std::vector<double> born::S_HCT;
//	std::vector<double> born::drOBC;
//	std::vector<double> born::U_ACE;
//	std::vector<double> born::W_ACE;
//	std::vector<double> born::S2_ACE;
//
//	surface::surface born::A;
//
//	double  born::offset = -0.09, born::probe = 1.4,
//		    born::electric = 332.0522173, born::off = 1.0E24,
//		    born::cut = 6.5E11;
//
//	double const  born::PI = 3.1415926535897932;
//
//	double born::A_solv = 0.0;
//
//	born::METHODS  born::METHOD = STILL;
//	born::SURFACES born::SURFACE = TINKER;
//	born::RADIUS   born::RSOLV = STD;
//
//	int born::DERIV = 0;
//
//	std::string born::method = "STILL";
//	std::string born::surface = "TINKER-Methode";
//	std::string born::rsolv = "TINKER standard vdW-Radien";
//	std::string born::RTYPE = "MACROMODEL";
//	bool born::IDEAL = true;
//	///////////////////////////////////////////////////////////////////
//}
//
//void GB::born::clean_methods()
//{
//	R_solv.clear();
//	V_solv.clear();
//	SA.clear();
//	S_HCT.clear();
//	drOBC.clear();
//	U_ACE.clear();
//	W_ACE.clear(),
//	S2_ACE.clear();
//
//	A.K.clear();
//}
//
//// Berechnung der Lösungmittel-zugänglichen Oberfläche nach der Methode im TINKER-Quelltext.
//// Die Lösungmittel-zugängliche Oberfläche der Atome wird nur grob angenähert.
////
//// Diese Methode benötigt die aktuellen Bornradien.
//void GB::born::get_SA_TINKER(void)
//{
//	SA.clear();
//	SA.resize(N);
//
//	for (size_t i = 0; i < alpha.size(); ++i)
//	{
//		if (alpha[i] != 0.0) SA[i] = 4.0*PI*pow(R_solv[i] + probe, 2.0)*pow(R_solv[i] / alpha[i], 6.0);
//	}
//}
//
//// get_SA_STILL() berechnet analytisch eine Näherung der Lösungsmittel-zugänglichen 
//// Oberflächen der Atome eines Moleküls.
////
////     Literatur :
////
////     W. C. Still, A. Tempczyk, R. C. Hawley, T. Hendrickson,
////     J. Am. Chem. Soc., 1990, 6127 - 6129.
////
//void GB::born::get_SA_STILL()
//{
//	std::vector<double>   p;
//	double   p12 = 0.95, p13 = 0.2, p14 = 0.35, Si, rij, bij,
//		     rik, ri, rj, pij, pijk;
//	size_t   oz, hn, indx;
//
//	SA.clear();
//	SA.resize(N);
//
//	p.resize(N);
//	for (size_t i = 0; i < N; ++i)
//	{
//		oz = OZ(i);
//
//		switch (oz)
//		{
//		case 1:
//			p[i] = 0.9822;
//			indx = bonds(i)[0];
//			if (OZ(indx) == 7)
//			{
//				hn = 0;
//				for (size_t j = 0; j < bonds(indx).size(); ++j)
//				{
//					if (OZ(bonds(indx)[j]) == 1) hn++;
//					if (OZ(bonds(indx)[j]) == 8)
//					{
//						p[i] = 0.9808;
//						break;
//					}
//				}
//				if (p[i] != 0.0)
//				{
//					p[i] = 1.2245;
//					if (hn == 3)
//						p[i] = 1.117;
//				}
//			}
//			break;
//		case 6:
//			p[i] = 2.39;
//
//			hn = 0;
//			for (size_t j = 0; j < bonds(i).size(); ++j)
//			{
//				if (OZ(bonds(i)[j]) == 1) hn++;
//			}
//			if (bonds(i).size() == 4)
//			{
//				if (hn == 3) p[i] = 0.9822;
//				if (hn == 2) p[i] = 1.0796;
//				if (hn == 1) p[i] = 1.3153;
//			}
//			else if (bonds(i).size() == 3)
//			{
//				p[i] = 1.6282;
//				for (size_t j = 0; j < bonds(i).size(); ++j)
//				{
//					if (OZ(bonds(i)[j]) == 8)
//					{
//						p[i] = 1.6526;
//						break;
//					}
//				}
//			}
//			else if (bonds(i).size() == 2)
//				p[i] = 1.4148;
//			break;
//		case 7:
//			p[i] = 1.297;
//			if (bonds(i).size() == 1)
//				p[i] = 1.0;
//			else if (bonds(i).size() == 2)
//			{
//				if (charge(i) >= 1.0)
//					p[i] = 1.65;
//				else
//				{
//					p[i] = 1.1;
//					for (size_t j = 0; j < bonds(i).size(); ++j)
//					{
//						if (OZ(bonds(i)[j]) == 1) p[i] = 1.1979;
//					}
//				}
//			}
//			else if (bonds(i).size() == 4) p[i] = 2.42;
//			break;
//		case 8:
//			p[i] = 1.0915;
//			if (bonds(i).size() == 1) p[i] = 1.0758;
//			break;
//		case 9:
//		case 17:
//		case 35:
//		case 53:
//			p[i] = 1.05;
//			break;
//		case 15:
//		case 16:
//			p[i] = 1.2915;
//			break;
//		default:
//			p[i] = 1.0;
//		}
//	}
//
//	double dAii_x, dAii_y, dAii_z, dAij_x, dAij_y, dAij_z;
//
//	for (size_t i = 0; i < N; ++i)
//	{
//		ri = R_solv[i];
//		Si = 4 * PI*pow(ri + probe, 2.0);
//
//		SA[i] = 1.0;
//
//		for (size_t j = 0; j < N; ++j)
//		{
//			rj = R_solv[j];
//			rij = (xyz(i) - xyz(j)).len();
//
//			if ((i != j) && (rij < rj + ri + 2 * probe))
//			{
//				bij = PI*(ri + probe)*(ri + rj + 2 * probe - rij)*(1.0 + (rj - ri) / rij);
//
//				pij = p14;
//				// überprüft, ob i mit j verbunden ist
//				for (size_t l = 0; l < bonds(i).size(); ++l)
//				{
//					if (j == bonds(i)[l])
//					{
//						pij = p12;
//						break;
//					}
//				}
//				// überprüft, ob i mit j über ein weiteres Atom verbunden ist
//				if (pij != p12)
//				{
//					for (size_t l = 0; l < bonds(i).size(); ++l)
//					{
//						for (size_t o = 0; o < bonds(bonds(i)[l]).size(); ++o)
//						{
//							if (j == bonds(bonds(i)[l])[o])
//							{
//								pij = p13;
//								break;
//							}
//						}
//						if (pij == p13) break;
//					}
//				}
//
//				pijk = 1.0;
//				for (size_t k = 0; k < bonds(i).size(); ++k)
//				{
//					if (j != bonds(i)[k])
//					{
//						rik = (xyz(i) - xyz(bonds(i)[k])).len();
//						pijk *= (rik*rik) / (rij*rij);
//					}
//				}
//				SA[i] *= (1.0 - (p[i] * pij*pijk*bij) / Si);
//
//				if (DERIV == 1)
//				{
//					dAii_x = -PI*(ri + probe)*(1.0 + ((ri + rj + 2 * probe)*(rj - ri)) / rij*rij)
//						/ ((Si / p[i] * pij) - bij);
//					dAii_y = dAii_x;
//					dAii_z = dAii_x;
//					dAii_x *= (xyz(i).x() - xyz(j).x()) / rij;
//					dAii_y *= (xyz(i).y() - xyz(j).y()) / rij;
//					dAii_z *= (xyz(i).z() - xyz(j).z()) / rij;
//					dAES[i].x() += dAii_x;
//					dAES[i].y() += dAii_y;
//					dAES[i].z() += dAii_z;
//					dAES[j].x() -= dAii_x;
//					dAES[j].y() -= dAii_y;
//					dAES[j].z() -= dAii_z;
//				}
//			}
//		}
//		SA[i] *= Si;
//		if (SA[i] < 0) SA[i] = 0.0;
//	}
//
//	if (DERIV == 1)
//	{
//		for (size_t i = 0; i < N; ++i)
//		{
//			dAES[i].x() *= SA[i];
//			dAES[i].y() *= SA[i];
//			dAES[i].z() *= SA[i];
//		}
//	}
//}
//
//// get_SA_GAUSS() berechnet die Lösungsmittel-zugänglichen Oberflächen der Atome
//// nach dem Theorem von Gauss-Bonnet.
//void GB::born::get_SA_GAUSS()
//{
//	surface::surface Atome;
//	std::vector<bool> omit;
//
//	Atome.K.resize(N);
//	omit.resize(N);
//	for (size_t i = 0; i < N; ++i)
//	{
//		Atome.K[i].m = xyz(i);
//		Atome.K[i].r = R_solv[i] + probe;
//
//		omit[i] = false;
//	}
//
//	SA.clear();
//	SA.resize(N);
//
//	for (size_t i = 0; i < N; ++i)
//	{
//		omit[i] = true;
//		if (i>0) omit[i - 1] = false;
//
//		SA[i] = Atome.GAUSS(Atome.K[i].m, Atome.K[i].r, omit, Atome.K, DERIV, dAES);
//	}
//}
//
//// GAUS() berechnet die Oberflächen von sich überschneidenden Kugeln nach dem
//// Theorem von Gauss-Bonnet.
//double surface::surface::GAUSS(scon::vect3d<double> M, double R, std::vector<bool> omit, std::vector<KUGEL> &K, int DERIV, coords::Representation_3D &dA)
//{
//	std::vector<size_t>   io;
//	std::vector<double>   dsq, b, gr, ther, bg, risq, ri;
//	std::vector<KUGEL>    J;
//
//	double   delta, eps, rmove, A;
//	double   arclen, exang;
//	size_t   jb, ib, maxarc(1000);
//	bool     moved;
//	int      ir(-1);
//
//	double   d, temp_d;
//	size_t   ink, temp_i;
//	scon::vect3d<double> temp_v;
//
//	if (R <= 0.0) return 0;
//
//	delta = 1.0E-8;
//	eps = 1.0E-8;
//	rmove = 1.0E-8;
//	moved = false;
//
//	A = 0.0;
//
//	jb = 0;
//	ib = 0;
//	arclen = 0.0;
//	exang = 0.0;
//
//
//	// Testet, welche der Kugeln in K sich mit Kugel {M,R} überschneidet
//	io.clear();
//	for (size_t i = 0; i < K.size(); ++i)
//	{
//		d = (M - K[i].m).len();
//		if ((std::max(R - K[i].r, K[i].r - R) <= d) && (d <= R + K[i].r))
//		if (!omit[i]) io.push_back(i);
//		if (omit[i]) ir = i;
//	}
//
//	// Paremeter der sich mit {M,R} überlappenden Kugeln werden bestimmt
//	J.clear();
//	gr.clear();
//	b.clear();
//	dsq.clear();
//
//	J.resize(io.size());
//	gr.resize(io.size());
//	b.resize(io.size());
//	dsq.resize(io.size());
//	for (size_t i = 0; i < io.size(); ++i)
//	{
//		// Verbindungsvektoren der Mittelpunkte
//		J[i].m = K[io[i]].m - M;
//		J[i].r = K[io[i]].r;
//
//		if (PERMUTATE)
//		{
//			temp_d = J[i].m.x();
//			J[i].m.setX(J[i].m.y());
//			J[i].m.setY(J[i].m.z());
//			J[i].m.setZ(temp_d);
//		}
//
//		dsq[i] = J[i].m.x() * J[i].m.x() + J[i].m.y() * J[i].m.y();
//
//		if ((dsq[i] < delta*delta) && ROUND)
//		{
//			dsq[i] = delta*delta;
//			J[i].m.setX(delta);
//			J[i].m.setY(0.0);
//		}
//		b[i] = J[i].m.len();
//		gr[i] = (b[i] * b[i] + (R + J[i].r)*(R - J[i].r)) / (2.0*R*b[i]);
//	}
//
//	// Keine der Kugeln in K schneidet sich it Kugel {M,R}
//	if (io.size() == 0) return 4.0 * PI * R*R;
//
//	// Geanu eine Kugel in K schneidet sich mit Kugel {M,R}
//	if (io.size() == 1)
//	{
//		A = (4.0*PI - 2.0 * PI * (1.0 - gr[0])) * R*R;
//
//		if (DERIV == 1)
//		{
//			double dt1 = 2.0*PI*R*R*(b[0] - R*R + pow(J[0].r, 2.0)) / (2.0*R*pow(b[0], 3.0));
//
//			if (!PERMUTATE)
//			{
//				dA[ir].x() -= J[0].m.x()*dt1;
//				dA[ir].y() -= J[0].m.y()*dt1;
//				dA[ir].z() -= J[0].m.z()*dt1;
//
//				dA[io[0]].x() += J[0].m.x()*dt1;
//				dA[io[0]].y() += J[0].m.y()*dt1;
//				dA[io[0]].z() += J[0].m.z()*dt1;
//			}
//			else
//			{
//				dA[ir].x() -= J[0].m.z()*dt1;
//				dA[ir].y() -= J[0].m.x()*dt1;
//				dA[ir].z() -= J[0].m.y()*dt1;
//
//				dA[io[0]].x() += J[0].m.z()*dt1;
//				dA[io[0]].y() += J[0].m.x()*dt1;
//				dA[io[0]].z() += J[0].m.y()*dt1;
//			}
//		}
//
//		return A;
//	}
//
//	// Anordnen der Kugeln nach ihrem Grad der Überlappung in aufsteigender Reihenfolge
//	ink = 0;
//	while (ink < gr.size() - 1)
//	{
//		ink += 1;
//		for (size_t i = 0; i < ink; ++i)
//		{
//			if (gr[i] > gr[ink])
//			{
//				temp_v = J[ink].m;
//				J[ink].m = J[i].m;
//				J[i].m = temp_v;
//
//				temp_d = J[ink].r;
//				J[ink].r = J[i].r;
//				J[i].r = temp_d;
//
//				temp_d = dsq[ink];
//				dsq[ink] = dsq[i];
//				dsq[i] = temp_d;
//
//				temp_d = b[ink];
//				b[ink] = b[i];
//				b[i] = temp_d;
//
//				temp_d = gr[ink];
//				gr[ink] = gr[i];
//				gr[i] = temp_d;
//
//				temp_i = io[ink];
//				io[ink] = io[i];
//				io[i] = temp_i;
//
//				ink = i;
//				break;
//			}
//		}
//	}
//
//	// Berechnet den Radius der Schnittkreise
//	ther.resize(gr.size());
//	bg.resize(gr.size());
//	risq.resize(gr.size());
//	ri.resize(gr.size());
//	for (size_t i = 0; i < gr.size(); ++i)
//	{
//		temp_d = gr[i] * R;
//		bg[i] = b[i] * temp_d;
//		risq[i] = R*R - temp_d*temp_d;
//		ri[i] = pow(risq[i], 0.5);
//		ther[i] = 0.5*PI - asin(std::min(1.0, std::max(-1.0, gr[i])));
//	}
//
//	// Auffinden der Grenzen der unzugänglichen Fläche auf Kugel {M,R}
//	double   cc, td;
//	bool     buried(false);
//	std::vector<bool>   omit2(io.size());
//
//	for (size_t k = 0; k < omit2.size(); ++k)
//		omit2[k] = false;
//	for (size_t k = 0; k < io.size() - 1; ++k)
//	{
//		if (!omit2[k])
//		{
//			for (size_t j = k + 1; j < io.size(); ++j)
//			{
//				if (!omit2[j])
//				{
//					cc = (J[k].m.x()*J[j].m.x() + J[k].m.y()*J[j].m.y() + J[k].m.z()*J[j].m.z()) / (b[k] * b[j]);
//					cc = acos(std::min(1.0, std::max(-1.0, cc)));
//					td = ther[k] + ther[j];
//					if (cc < td)
//					{
//						if (cc + ther[j] < ther[k]) omit2[j] = true;
//						else if (cc > delta) if (2.0*PI - cc <= td) buried = true;
//					}
//				}
//			}
//		}
//		if (buried) break;
//	}
//
//	if (!buried)
//	{
//		std::vector<double>   arcf, arci, ex, ux, uy, uz;
//		std::vector<size_t>   lt, kent, kout, ider;
//		std::vector<int>   sign_yder;
//
//		double   t1, dk, cosine, dsqj, tb, txb, tyb,
//			tr, tr2, txr, tyr, tk1, tk2, thec, the,
//			ti, tf, arcsum, t;
//		double   axx, axy, axz, ayx, ayy, azx, azy, azz;
//		double   uxj, uyj, uzj;
//		size_t   narc, ni, nj;
//		bool     top;
//
//		kent.clear();
//		kout.clear();
//
//		for (size_t k = 0; k < io.size(); ++k)
//		{
//			if (!omit2[k])
//			{
//				omit2[k] = true;
//				narc = 0;
//				dk = pow(dsq[k], 0.5);
//				top = false;
//
//				t1 = J[k].m.z() / (b[k] * dk);
//				axx = J[k].m.x() * t1;
//				axy = J[k].m.y() * t1;
//				axz = dk / b[k];
//				ayx = J[k].m.y() / dk;
//				ayy = J[k].m.x() / dk;
//				azx = J[k].m.x() / b[k];
//				azy = J[k].m.y() / b[k];
//				azz = J[k].m.z() / b[k];
//
//				arcf.clear();
//				arci.clear();
//				ex.clear();
//				lt.clear();
//
//				if (DERIV == 1)
//				{
//					ux.clear();
//					uy.clear();
//					uz.clear();
//					ider.clear();
//					sign_yder.clear();
//
//					ider.resize(io.size());
//					sign_yder.resize(io.size());
//					ux.resize(io.size());
//					uy.resize(io.size());
//					uz.resize(io.size());
//				}
//
//				for (size_t j = 0; j < io.size(); ++j)
//				{
//					if (!omit2[j])
//					{
//						uxj = J[j].m.x()*axx + J[j].m.y()*axy - J[j].m.z()*axz;
//						uyj = J[j].m.y()*ayy - J[j].m.x()*ayx;
//						uzj = J[j].m.x()*azx + J[j].m.y()*azy + J[j].m.z()*azz;
//						cosine = std::min(1.0, std::max(-1.0, uzj / b[j]));
//
//						if (acos(cosine) < (ther[k] + ther[j]))
//						{
//							dsqj = uxj*uxj + uyj*uyj;
//							tb = uzj*gr[k] * R - bg[j];
//							txb = uxj * tb;
//							tyb = uyj * tb;
//							td = ri[k] * dsqj;
//							tr2 = risq[k] * dsqj - tb*tb;;
//							tr2 = std::max(eps, tr2);
//							tr = pow(tr2, 0.5);
//							txr = uxj * tr;
//							tyr = uyj * tr;
//
//							tb = (txb + tyr) / td;
//							tb = std::min(1.0, std::max(-1.0, tb));
//							tk1 = acos(tb);
//							if (tyb - txr < 0.0) tk1 = 2.0*PI - tk1;
//							tb = (txb - tyr) / td;
//							tb = std::min(1.0, std::max(-1.0, tb));
//							tk2 = acos(tb);
//
//							if (tyb + txr < 0.0) tk2 = 2.0*PI - tk2;
//							thec = (R*R*uzj - gr[k] * R*bg[j]) / (ri[k] * ri[j] * b[j]);
//							if (std::max(thec, -1.0*thec) < 1.0) the = -1.0 * acos(thec);
//							else if (thec >= 1.0) the = 0.0;
//							else if (thec <= -1.0) the = -1.0*PI;
//							cosine = std::min(1.0, std::max(-1.0, (uzj*gr[k] * R - uxj*ri[k]) / (b[j] * R)));
//
//							if (((acos(cosine) - ther[j])*(tk2 - tk1)) <= 0.0)
//							{
//								ti = tk2;
//								tf = tk1;
//							}
//							else
//							{
//								// Fehler in TINKER-Routine surfatom (Born-Radien Berechnung): else-Anweisung wie if-Anweisung
//								ti = tk2;
//								tf = tk1;
//								// in TINKER-Routine surface1, die zur Berechnung der Ableitungen benutzt wird,
//								// steht die korrekte Anweisung
//								/*if (DERIV == 1)
//								{
//								ti = tk1;
//								tf = tk2;
//								}*/
//							}
//							narc += 1;
//							if (tf <= ti)
//							{
//								arcf.push_back(tf);
//								arci.push_back(0.0);
//								tf = 2.0*PI;
//								lt.push_back(j);
//								ex.push_back(the);
//								top = true;
//								narc += 1;
//							}
//							arcf.push_back(tf);
//							arci.push_back(ti);
//							lt.push_back(j);
//							ex.push_back(the);
//
//							if (DERIV == 1)
//							{
//								ux[j] = uxj;
//								uy[j] = uyj;
//								uz[j] = uzj;
//							}
//						}
//					}
//				}
//				omit2[k] = false;
//
//				if (narc <= 0)
//				{
//					arcsum = 2.0*PI;
//					ib = ib + 1;
//				}
//				else
//				{
//					ink = 0;
//					while (ink < arci.size() - 1)
//					{
//						ink += 1;
//						for (size_t i = 0; i < ink; ++i)
//						{
//							if (arci[i]>arci[ink])
//							{
//								temp_d = arci[ink];
//								arci[ink] = arci[i];
//								arci[i] = temp_d;
//
//								temp_d = arcf[ink];
//								arcf[ink] = arcf[i];
//								arcf[i] = temp_d;
//
//								temp_d = ex[ink];
//								ex[ink] = ex[i];
//								ex[i] = temp_d;
//
//								temp_i = lt[ink];
//								lt[ink] = lt[i];
//								lt[i] = temp_i;
//
//								ink = i;
//								break;
//							}
//						}
//					}
//
//					arcsum = arci[0];
//					t = arcf[0];
//					ni = 0;
//
//					if (narc > 1)
//					{
//						for (size_t j = 1; j < narc; ++j)
//						{
//							if (t < arci[j])
//							{
//								arcsum = arcsum + arci[j] - t;
//								exang = exang + ex[ni];
//								jb = jb + 1;
//								temp_i = lt[ni];
//								if (DERIV == 1) ider[temp_i] += 1;
//								if (DERIV == 1) sign_yder[temp_i] += 1;
//								kent.push_back(maxarc*temp_i + k);
//								temp_i = lt[j];
//								if (DERIV == 1) ider[temp_i] += 1;
//								if (DERIV == 1) sign_yder[temp_i] -= 1;
//								kout.push_back(maxarc*k + temp_i);
//							}
//							if (arcf[j] >= t)
//							{
//								t = arcf[j];
//								ni = j;
//							}
//						}
//					}
//					arcsum = arcsum + 2.0*PI - t;
//
//					if (!top)
//					{
//						exang = exang + ex[ni];
//						jb = jb + 1;
//						temp_i = lt[ni];
//						if (DERIV == 1) ider[temp_i] += 1;
//						if (DERIV == 1) sign_yder[temp_i] += 1;
//						kent.push_back(maxarc*temp_i + k);
//						temp_i = lt[0];
//						if (DERIV == 1) ider[temp_i] += 1;
//						if (DERIV == 1) sign_yder[temp_i] -= 1;
//						kout.push_back(maxarc*k + temp_i);
//					}
//
//					// Ableitungen der Oberflächen werden bestimmt
//					if (DERIV == 1)
//					{
//						double   rcn, uzl, gl, bgl, bsql, risql, wxlsq, wxl, p, v, tl,
//							tm, gk, deal, decl, dtkal, dtkcl, s, dtlal, dtlcl,
//							gaca, gacb, faca, facb, facc, dax, day, daz;
//
//						for (size_t l = 0; l < io.size(); ++l)
//						{
//							if (ider[l] != 0)
//							{
//								rcn = ider[l] * R*R;
//								ider[l] = 0;
//								uzl = uz[l];
//								gl = gr[l] * R;
//								gk = gr[k] * R;
//								bgl = bg[l];
//								bsql = pow(b[l], 2.0);
//								risql = risq[l];
//								wxlsq = bsql - uzl*uzl;
//								wxl = pow(std::fabs(wxlsq), 0.5);
//								p = bgl - gk * uzl;
//								v = risq[k] * wxlsq - p*p;
//								if (v < 0.0) v *= -1.0;
//								v = std::max(eps, v);
//								v = pow(v, 0.5);
//								tl = R*(gk*(bgl - bsql) + uzl*(bgl - R*R)) / (v*risql*bsql);
//								deal = -wxl*tl;
//								decl = -uzl*tl - R / v;
//								dtkal = (wxlsq - p) / (wxl*v);
//								dtkcl = (uzl - gk) / v;
//								s = gk*b[l] - gl * uzl;
//								tl = 2.0*gk - uzl;
//								tm = R*R - bgl;
//								dtlal = -(risql*wxlsq*b[l] * tl - s*(wxlsq*tm + risql*bsql)) / (risql*wxl*bsql*v);
//								dtlcl = -(risql*b[l] * (uzl*tl - bgl) - uzl*tm*s) / (risql*bsql*v);
//								gaca = rcn*(deal - (gk*dtkal - gl*dtlal) / R) / wxl;
//								gacb = (gk - uzl*gl / b[l])*sign_yder[l] * R / wxlsq;
//								sign_yder[l] = 0;
//
//								faca = ux[l] * gaca - uy[l] * gacb;
//								facb = uy[l] * gaca + ux[l] * gacb;
//								facc = rcn*(decl - (gk*dtkcl - gl*dtlcl) / R);
//								dax = axx*faca - ayx*facb + azx*facc;
//								day = axy*faca + ayy*facb + azy*facc;
//								daz = azz*facc - axz*faca;
//
//								if (!PERMUTATE)
//								{
//									dA[ir].x() += dax;
//									dA[ir].y() += day;
//									dA[ir].z() += daz;
//
//									dA[io[l]].x() -= dax;
//									dA[io[l]].y() -= day;
//									dA[io[l]].z() -= daz;
//								}
//								else
//								{
//									dA[ir].x() += daz;
//									dA[ir].y() += dax;
//									dA[ir].z() += day;
//
//									dA[io[l]].x() -= daz;
//									dA[io[l]].y() -= dax;
//									dA[io[l]].z() -= day;
//								}
//							}
//						}
//					}
//				}
//
//				arclen = arclen + gr[k] * arcsum;
//
//				if (DERIV == 1)
//				{
//					double dt1 = arcsum*R*R*(b[k] * b[k] - R*R + pow(J[k].r, 2.0)) / (2.0*R*pow(b[k], 3.0));
//
//					if (!PERMUTATE)
//					{
//						dA[ir].x() -= J[k].m.x()*dt1;
//						dA[ir].y() -= J[k].m.y()*dt1;
//						dA[ir].z() -= J[k].m.z()*dt1;
//
//						dA[io[k]].x() += J[k].m.x()*dt1;
//						dA[io[k]].y() += J[k].m.y()*dt1;
//						dA[io[k]].z() += J[k].m.z()*dt1;
//					}
//					else
//					{
//						dA[ir].x() -= J[k].m.z()*dt1;
//						dA[ir].y() -= J[k].m.x()*dt1;
//						dA[ir].z() -= J[k].m.y()*dt1;
//
//						dA[io[k]].x() += J[k].m.z()*dt1;
//						dA[io[k]].y() += J[k].m.x()*dt1;
//						dA[io[k]].z() += J[k].m.y()*dt1;
//					}
//				}
//			}
//		}
//		// Findet die Anzahl der unabhängigen Bogensegmente und überprüft, ob sie verbunden sind
//		nj = 0;
//		bool   repeat, jump150(false);
//		size_t   m;
//
//		for (size_t k = 0; k < jb; ++k)
//		{
//			if (kout[k] != 0)
//			{
//				repeat = true;
//				temp_i = k;
//
//				while (repeat)
//				{
//					repeat = false;
//
//					m = kout[temp_i];
//					kout[temp_i] = 0;
//					nj = nj + 1;
//
//					for (size_t ii = 0; ii < jb; ++ii)
//					{
//						if (m == kent[ii])
//						{
//							if (ii == k)
//							{
//								ib = ib + 1;
//								break;
//							}
//							else
//							{
//								temp_i = ii;
//								repeat = true;
//								break;
//							}
//
//						}
//					}
//				}
//			}
//		}
//
//		A = ib*2.0*PI + exang + arclen;
//		A = R*R * (A - std::floor(A / (4.0*PI))*4.0*PI);
//	}
//	return A;
//};