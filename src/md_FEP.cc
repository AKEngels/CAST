#include "md.h"

// If FEP calculation is requested: calculate lambda values for each window
// and print the scaling factors for van-der-Waals and electrostatics for each window
void md::simulation::fepinit(void)
{
  if (Config::get().general.verbosity > 1)   // print warning if vdw-softcore potential is shifted so strongly that is has no minimum
  {
    if (Config::get().general.energy_interface == config::interface_types::FORCEFIELD)
    {
      if (Config::get().fep.ljshift > 1.0 / ((1 - Config::get().fep.dlambda) * (1 - Config::get().fep.dlambda)))
      {
        std::cout << "WARNING! You should choose a smaller value for vshift!\n\n";
      }
    }
  }
  init();
  // center and temp var
  coordobj.move_all_by(-coordobj.center_of_geometry());
  double  increment;
  Config::set().md.fep = true;
  FEPsum = 0.0;
  FEPsum_back = 0.0;
  FEPsum_SOS = 0.0;
  FEPsum_BAR = 0.0;
  // fill vector with scaling increments
  increment = Config::get().fep.lambda / Config::get().fep.dlambda;
  //std::cout << Config::get().fep.lambda << "   " << Config::get().fep.dlambda << std::endl;
  std::cout << "Number of FEP windows:  " << increment << std::endl;
  coordobj.getFep().window.resize(std::size_t(increment) + 1); //because of backward transformation one window more necessary
  coordobj.getFep().window[0].step = 0;

  // calculate all lambda values for every window
  for (auto i = 0u; i < coordobj.getFep().window.size(); i++) {
    double lambda = i * Config::get().fep.dlambda;  // lambda
    if (lambda < Config::get().fep.eleccouple) coordobj.getFep().window[i].ein = 0;
    else coordobj.getFep().window[i].ein = 1 - (1 - lambda) / (1 - Config::get().fep.eleccouple);
    if (lambda > Config::get().fep.vdwcouple) coordobj.getFep().window[i].vin = 1;
    else coordobj.getFep().window[i].vin = lambda / Config::get().fep.vdwcouple;
    if (lambda > 1 - Config::get().fep.eleccouple) coordobj.getFep().window[i].eout = 0;
    else coordobj.getFep().window[i].eout = 1 - (lambda) / (1 - Config::get().fep.eleccouple);
    if (lambda < 1 - Config::get().fep.vdwcouple) coordobj.getFep().window[i].vout = 1;
    else coordobj.getFep().window[i].vout = (1 - lambda) / Config::get().fep.vdwcouple;

    double dlambda = (i + 1) * Config::get().fep.dlambda;  // lambda + dlambda

    if (dlambda < Config::get().fep.eleccouple) coordobj.getFep().window[i].dein = 0;
    else if (dlambda > 1) coordobj.getFep().window[i].dein = 1;
    else coordobj.getFep().window[i].dein = 1 - (1 - dlambda) / (1 - Config::get().fep.eleccouple);
    if (dlambda > Config::get().fep.vdwcouple) coordobj.getFep().window[i].dvin = 1;
    else coordobj.getFep().window[i].dvin = dlambda / Config::get().fep.vdwcouple;
    if (dlambda > 1 - Config::get().fep.eleccouple) coordobj.getFep().window[i].deout = 0;
    else coordobj.getFep().window[i].deout = 1 - (dlambda) / (1 - Config::get().fep.eleccouple);
    if (dlambda < 1 - Config::get().fep.vdwcouple) coordobj.getFep().window[i].dvout = 1;
    else if (dlambda > 1) coordobj.getFep().window[i].dvout = 0;
    else coordobj.getFep().window[i].dvout = (1 - dlambda) / Config::get().fep.vdwcouple;

    double mlambda = (int)(i - 1) * Config::get().fep.dlambda;  // lambda - dlambda
    if (mlambda < Config::get().fep.eleccouple) coordobj.getFep().window[i].mein = 0;
    else coordobj.getFep().window[i].mein = 1 - (1 - mlambda) / (1 - Config::get().fep.eleccouple);
    if (mlambda > Config::get().fep.vdwcouple) coordobj.getFep().window[i].mvin = 1;
    else if (mlambda < 0) coordobj.getFep().window[i].mvin = 0;
    else coordobj.getFep().window[i].mvin = mlambda / Config::get().fep.vdwcouple;
    if (mlambda > 1 - Config::get().fep.eleccouple) coordobj.getFep().window[i].meout = 0;
    else if (mlambda < 0) coordobj.getFep().window[i].meout = 1;
    else coordobj.getFep().window[i].meout = 1 - (mlambda) / (1 - Config::get().fep.eleccouple);
    if (mlambda < 1 - Config::get().fep.vdwcouple) coordobj.getFep().window[i].mvout = 1;
    else coordobj.getFep().window[i].mvout = (1 - mlambda) / Config::get().fep.vdwcouple;

  }// end of loop

  // clear FEP output vector and print lambvda values
  std::ofstream fepclear;
  fepclear.open("alchemical.txt");
  fepclear.close();
  coordobj.getFep().window[0].step = 0;
  std::cout << "FEP Coupling Parameters:" << std::endl;
  std::cout << std::endl;
  std::cout << "Van-der-Waals Coupling: " << std::endl;
  for (std::size_t i = 0; i < coordobj.getFep().window.size(); i++)
  {
    std::cout << std::setprecision(4) << std::setw(8) << coordobj.getFep().window[i].mvout << std::setw(8) << coordobj.getFep().window[i].vout << std::setw(8) << coordobj.getFep().window[i].dvout << std::setw(8) << coordobj.getFep().window[i].mvin << std::setw(8) << coordobj.getFep().window[i].vin << std::setw(8) << coordobj.getFep().window[i].dvin << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Electrostatic Coupling:" << std::endl;
  for (std::size_t i = 0; i < coordobj.getFep().window.size(); i++)
  {
    std::cout << std::setprecision(4) << std::setw(8) << coordobj.getFep().window[i].meout << std::setw(8) << coordobj.getFep().window[i].eout << std::setw(8) << coordobj.getFep().window[i].deout << std::setw(8) << coordobj.getFep().window[i].mein << std::setw(8) << coordobj.getFep().window[i].ein << std::setw(8) << coordobj.getFep().window[i].dein << std::endl;
  }

}

// Calculation of ensemble average and free energy change for each step if FEP calculation is performed
// calculation can be improved if at every step the current averages are stored
// currently calculation is performed at the end of each window
void md::simulation::freecalc()
{
  std::size_t iterator(0U), k(0U);
  // set conversion factors (conv) and constants (boltzmann, avogadro)
  double de_ensemble, de_ensemble_back, de_ensemble_half = 0, de_ensemble_back_half = 0, temp_avg = 0.0, boltz = 1.3806488E-23, avogad = 6.022E23, conv = 4184.0;
  // calculate ensemble average for current window
  for (std::size_t i = 0; i < coordobj.getFep().fepdata.size(); i++)
  {                // for every conformation in window
    iterator += 1;
    k = 0;
    de_ensemble = de_ensemble_back = temp_avg = 0.0;
    double exponent = -1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * coordobj.getFep().fepdata[i].dE / avogad;
    double exponent_back = 1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * coordobj.getFep().fepdata[i].dE_back / avogad;
    coordobj.getFep().fepdata[i].de_ens = exp(exponent);
    coordobj.getFep().fepdata[i].de_ens_back = exp(exponent_back);
    double de_ens_half = exp(exponent / 2);
    double de_ens_back_half = exp(exponent_back / 2);
    for (k = 0; k <= i; k++) {
      temp_avg += coordobj.getFep().fepdata[k].T;
      de_ensemble += coordobj.getFep().fepdata[k].de_ens;
      de_ensemble_back += coordobj.getFep().fepdata[k].de_ens_back;
    }
    de_ensemble = de_ensemble / iterator;
    de_ensemble_back = de_ensemble_back / iterator;
    de_ensemble_half += de_ens_half;
    de_ensemble_back_half += de_ens_back_half;
    temp_avg = temp_avg / iterator;
    coordobj.getFep().fepdata[i].dG = -1 * std::log(de_ensemble) * temp_avg * boltz * avogad / conv;
    coordobj.getFep().fepdata[i].dG_back = std::log(de_ensemble_back) * temp_avg * boltz * avogad / conv;
  }// end of main loop

  de_ensemble_half = de_ensemble_half / coordobj.getFep().fepdata.size();
  de_ensemble_back_half = de_ensemble_back_half / coordobj.getFep().fepdata.size();
  if (de_ensemble_back_half == 1)  dG_SOS = 0;
  else dG_SOS = -1 * std::log(de_ensemble_v_SOS / de_ensemble_back_half) * temp_avg * boltz * avogad / conv;
  de_ensemble_v_SOS = de_ensemble_half; //de_ensemble_half is needed for SOS-calculation in next step

  // calculate final free energy change for the current window
  this->FEPsum += coordobj.getFep().fepdata[coordobj.getFep().fepdata.size() - 1].dG;
  this->FEPsum_back += coordobj.getFep().fepdata[coordobj.getFep().fepdata.size() - 1].dG_back;
  this->FEPsum_SOS += dG_SOS;
}

void md::simulation::bar(int current_window)
{
  double boltz = 1.3806488E-23, avogad = 6.022E23, conv = 4184.0;
  double w, w_back;  // weighting function
  double dG_BAR = dG_SOS;  // start value for iteration
  double c; // constant C

  if (Config::get().general.verbosity > 3)
  {
    std::cout << "Start solution of BAR equation from dG_SOS: " << dG_SOS << "\n";
  }

  auto count_iterations{ 0u };
  do    // iterative solution for BAR equation
  {
    c = dG_BAR;
    double ensemble = 0;
    double ensemble_back = 0;
    double temp_avg = 0;  // average temperature
    for (std::size_t i = 0; i < coordobj.getFep().fepdata.size(); i++) // for every conformation in window
    {
      w = 2 / (exp(1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * ((coordobj.getFep().fepdata[i].dE - c) / 2) / avogad) + exp(-1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * ((coordobj.getFep().fepdata[i].dE - c) / 2) / avogad));
      double ens = w * exp(-1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * (coordobj.getFep().fepdata[i].dE / 2) / avogad);
      w_back = 2 / (exp(1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * ((coordobj.getFep().fepdata[i].dE_back - c) / 2) / avogad) + exp(-1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * ((coordobj.getFep().fepdata[i].dE_back - c) / 2) / avogad));
      double ens_back = w_back * exp(1 / (boltz * coordobj.getFep().fepdata[i].T) * conv * (coordobj.getFep().fepdata[i].dE_back / 2) / avogad);
      ensemble += ens;
      ensemble_back += ens_back;
      temp_avg += coordobj.getFep().fepdata[i].T;
    }
    ensemble = ensemble / coordobj.getFep().fepdata.size();       // calculate averages
    ensemble_back = ensemble_back / coordobj.getFep().fepdata.size();
    temp_avg = temp_avg / coordobj.getFep().fepdata.size();

    if (current_window == 0)  dG_BAR = 0;                            // calculate dG for current window
    else dG_BAR = -1 * std::log(de_ensemble_v_BAR / ensemble_back) * temp_avg * boltz * avogad / conv;
    if (Config::get().general.verbosity > 3)
    {
      std::cout << "dG_BAR: " << dG_BAR << "\n";
    }
    de_ensemble_v_BAR = ensemble;  // this is needed in next step

    count_iterations += 1;
  } while (fabs(c - dG_BAR) > 0.001 && count_iterations < 5000);  // 0.001 = convergence threshold and 5000 = maximum number of iterations (maybe later define by user?)
  this->FEPsum_BAR += dG_BAR;
}

// write the output FEP calculations
void md::simulation::freewrite(int i)
{
  std::ofstream fep("alchemical.txt", std::ios_base::app);
  std::ofstream res("FEP_Results.txt", std::ios_base::app);

  // standard forward output
  if (i * Config::get().fep.dlambda == 0 && this->prod == false)
  {
    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << "0" << std::setw(10) << "0";
  }
  // equilibration is performed
  if (this->prod == false) {
    fep << "Equilibration for Lambda =  " << i * Config::get().fep.dlambda <<
      "  and dLambda =  " << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::endl;
  }
  // production run is performed
  else if (this->prod == true) {
    fep << "Starting new data collection with values:  " <<
      i * Config::get().fep.dlambda << "   " << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::endl;
  }
  // write output to alchemical.txt
  for (std::size_t k = 0; k < coordobj.getFep().fepdata.size(); k++) {
    if (k % Config::get().fep.freq == 0) {
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_c_l0;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_c_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_c_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_vdw_l0;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_vdw_l1;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].e_vdw_l2;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].T;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].dE;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].dG;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].dE_back;
      fep << std::fixed << std::right << std::setprecision(4) << std::setw(15) << coordobj.getFep().fepdata[k].dG_back;
      fep << std::endl;
    }
    if (Config::get().general.verbosity > 3u)
    {
      std::cout << "Coulomb: " << coordobj.getFep().fepdata[k].e_c_l2 - coordobj.getFep().fepdata[k].e_c_l1 << ", vdW: " << coordobj.getFep().fepdata[k].e_vdw_l2 - coordobj.getFep().fepdata[k].e_vdw_l1 << "\n";
    }
  }

  // at the end of production data in alchemical.txt sum up the results and print the before the new window starts
  if (this->prod == true) {
    fep << "Free energy change for the current window:  ";
    fep << coordobj.getFep().fepdata[coordobj.getFep().fepdata.size() - 1].dG << std::endl;
    fep << "Total free energy change until current window:  " << FEPsum << std::endl;

    res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << FEPsum_back << std::right << std::setprecision(4) << std::setw(10) << FEPsum_SOS << std::right << std::setprecision(4) << std::setw(10) << FEPsum_BAR << std::endl;
    double rounded = std::stod(std::to_string(i * Config::get().fep.dlambda)); // round necessary when having a number of windows that is can't be expressed exactly in decimal numbers
    if (rounded < 1) {
      res << std::fixed << std::right << std::setprecision(4) << std::setw(10) << (i * Config::get().fep.dlambda) + Config::get().fep.dlambda << std::setw(10) << FEPsum;
    }
  }
}
#ifdef USE_PYTHON
std::string md::simulation::get_pythonpath()
{
  std::string pythonpaths_str = Py_GetPath();
  std::string path;
#ifdef __unix__
  std::vector<std::string> pythonpaths = split(pythonpaths_str, ':');
#elif defined(_WIN32) || defined(WIN32)
  std::vector<std::string> pythonpaths = split(pythonpaths_str, ';');
#endif
  path = "import sys\n";
  for (auto p : pythonpaths)  //keep pythonpath of system
  {
    path += "sys.path.append('" + p + "')\n";
  }
  path += "sys.path.append('" + get_python_modulepath("matplotlib") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("fractions") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("csv") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("atexit") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("calendar") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("Tkinter") + "')\n";
  path += "sys.path.append('" + get_python_modulepath("FileDialog") + "')\n";
  return path;
}

std::vector<double> md::simulation::fepanalyze(std::vector<double> dE_pots, int window)
{
  std::string add_path = get_pythonpath();

  PyObject* modul, * funk, * prm, * ret, * pValue;

  // create python list with dE_pot_bac
  PyObject* E_pot_backs = PyList_New(coordobj.getFep().fepdata.size());
  for (std::size_t k = 0; k < coordobj.getFep().fepdata.size(); k++) {
    pValue = PyFloat_FromDouble(coordobj.getFep().fepdata[k].dE_back);
    PyList_SetItem(E_pot_backs, k, pValue);
  }

  if (window > 0)   // no output for 0th window
  {
    // create python list with dE_pot from last run
    PyObject* E_pots = PyList_New(coordobj.getFep().fepdata.size());
    for (std::size_t k = 0; k < coordobj.getFep().fepdata.size(); k++) {
      pValue = PyFloat_FromDouble(dE_pots[k]);
      PyList_SetItem(E_pots, k, pValue);
    }

    PySys_SetPath((char*)"./python_modules"); //set path
    const char* c = add_path.c_str();  //add paths pythonpath
    PyRun_SimpleString(c);

    modul = PyImport_ImportModule("FEP_analysis"); //import module 
    if (modul)
    {
      funk = PyObject_GetAttrString(modul, "plot_histograms_and_calculate_overlap"); //create function
      prm = Py_BuildValue("OOi", E_pots, E_pot_backs, window); //give parameters
      ret = PyObject_CallObject(funk, prm);  //call function with parameters
      std::string result_str = PyString_AsString(ret); //convert result to a C++ string
      if (result_str == "error")
      {
        std::cout << "An error occured during running python module 'FEP_analysis'\n";
      }
      else  // python function was successfull
      {
        float result = std::stof(result_str);  // convert result to float
        std::ofstream overlap("overlap.txt", std::ios_base::app);
        overlap << "Window " << window << ": " << result * 100 << " %\n";
      }
    }
    else
    {
      throw std::runtime_error("Error: module 'FEP_analysis' not found!");
    }
    //delete PyObjects
    Py_DECREF(prm);
    Py_DECREF(ret);
    Py_DECREF(funk);
    Py_DECREF(modul);
    Py_DECREF(pValue);
    Py_DECREF(E_pots);
    Py_DECREF(E_pot_backs);
  }

  dE_pots.clear();  // save dE_pot for next run and return them
  for (std::size_t k = 0; k < coordobj.getFep().fepdata.size(); k++) {
    dE_pots.push_back(coordobj.getFep().fepdata[k].dE);
  }
  return dE_pots;
}
#endif

// perform FEP calculation if requested
void md::simulation::feprun()
{
  if (Config::get().fep.analyze)
  {
    std::remove("overlap.txt");
  }
  std::vector<double> dE_pots;

  for (auto i(0U); i < coordobj.getFep().window.size(); ++i)  //for every window
  {
    std::cout << "Lambda:  " << i * Config::get().fep.dlambda << "\n";
    coordobj.getFep().window[0U].step = static_cast<int>(i);
    coordobj.getFep().fepdata.clear();
    // equilibration run for window i
    Config::set().md.num_steps = Config::get().fep.equil;
    integrate(true);
    // write output for equlibration and clear fep vector
    this->prod = false;
    freewrite(i);
    coordobj.getFep().fepdata.clear();
    // production run for window i
    Config::set().md.num_steps = Config::get().fep.steps;
    integrate(true);
    this->prod = true;
    // calculate free energy change for window and write output
    freecalc();
    if (Config::get().fep.bar == true) bar(i);
    freewrite(i);

    if (Config::get().fep.analyze)
    {
#ifdef USE_PYTHON
      dE_pots = fepanalyze(dE_pots, i);
#else
      std::cout << "Analyzing is not possible without python!\n";
#endif
    }
  }
}// end of main window loop
