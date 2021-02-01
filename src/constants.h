// Constants
#pragma once

namespace constants
{
// These constants define all units in CAST.
  constexpr double N_avogadro = 6.02214076e23; // 1/mol, Avogadros constant
  constexpr double pi = 3.1415926535897932384626433832795029;
  constexpr double e = 2.71828182845904523536028747135266249775724709369995;
  constexpr double gamma = 0.5772156649015328606065120900824024L; // Eulers constant
  constexpr double ln_2 = 0.6931471805599453094172321214581766; // the natural logarithm of 2 in high precision
  constexpr double h_bar_SI_units = 1.054571817*1e-34; // in J/K
  constexpr double h_bar_gaussian_units = 1. / (2. * pi) * 4.135667662e-15; // in eV * s
  constexpr double epsilon_0 = 8.854187e-12; // in units_ farad per meter (F/m)
  constexpr double boltzmann_constant_kb_gaussian_units = 8.617333262145e-5; //  in gauss einheiten // Dustin July19: is in eV/K
  constexpr double eV2kcal_mol(23.06035);
  constexpr double gas_constant_R_CASTunits = 0.83144725; // This is R in g*angstrom*angstrom/(picosecond*picosecond*Kelvin*mol)
  constexpr double gas_constant_R_kcal_per_mol_kelvin = 1.9872066e-3; // [kcal/(K*mol)]
  constexpr double au2kcal_mol(627.5095);
  constexpr double cal2joules(4.184);
  constexpr double speed_of_light_m_per_s(299792458); // meters per second
  constexpr double speed_of_light_cm_per_s(speed_of_light_m_per_s*100);
  //
  constexpr double kcal_mol2ev(1.0 / eV2kcal_mol);
  constexpr double kcal_mol2au(1.0 / au2kcal_mol);
  constexpr double joules2cal(1.0/ cal2joules);
  constexpr double boltzmann_constant_kb_SI_units = boltzmann_constant_kb_gaussian_units * eV2kcal_mol * cal2joules * 1000. / N_avogadro; //J/K
}

namespace mathFunctions
{


  /////////////////
  // Auxiliary Functions
  /////////////////

  // Gives DiGamma Function value, via Stackoverflow
  inline long double digammal(long double x)
  {
    using namespace constants;
    /* force into the interval 1..3 */
    if (x < 0.0L)
      return digammal(1.0L - x) + pi / tanl(pi * (1.0L - x));	/* reflection formula */
    else if (x < 1.0L)
      return digammal(1.0L + x) - 1.0L / x;
    else if (x == 1.0L)
      return -::constants::gamma;
    else if (x == 2.0L)
      return 1.0L - ::constants::gamma;
    else if (x == 3.0L)
      return 1.5L - ::constants::gamma;
    else if (x > 3.0L)
      /* duplication formula */
      return 0.5L * (digammal(x / 2.0L) + digammal((x + 1.0L) / 2.0L)) + ln_2;
    else
    {
      /* Just for your information, the following lines contain
      * the Maple source code to re-generate the table that is
      * eventually becoming the Kncoe[] array below
      * interface(prettyprint=0) :
      * Digits := 63 :
      * r := 0 :
      *
      * for l from 1 to 60 do
      * 	d := binomial(-1/2,l) :
      * 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
      * 	evalf(r) ;
      * 	print(%,evalf(1+Psi(1)-r)) ;
      *o d :
      *
      * for N from 1 to 28 do
      * 	r := 0 :
      * 	n := N-1 :
      *
      *	for l from iquo(n+3,2) to 70 do
      *		d := 0 :
      *		for s from 0 to n+1 do
      *		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
      *		od :
      *		if 2*l-n > 1 then
      *		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
      *		fi :
      *	od :
      *	print(evalf((-1)^n*2*r)) ;
      *od :
      *quit :
      */
      static long double Kncoe[] = { .30459198558715155634315638246624251L,
        .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
        .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
        .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
        .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
        .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
        .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
        .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
        .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
        .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
        .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
        .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
        .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
        .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
        .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L };

      long double Tn_1 = 1.0L;	/* T_{n-1}(x), started at n=1 */
      long double Tn = x - 2.0L;	/* T_{n}(x) , started at n=1 */
      long double resul = Kncoe[0] + Kncoe[1] * Tn;

      x -= 2.0L;

      for (auto n = 2u; n < sizeof(Kncoe) / sizeof(long double); n++)
      {
        const long double Tn1 = 2.0L * x * Tn - Tn_1;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
        resul += Kncoe[n] * Tn1;
        Tn_1 = Tn;
        Tn = Tn1;
      }
      return resul;
    }
  }


}
