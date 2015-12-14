/*

 Lit.
 [1] Kaestner, Sherwood, J. Chem. Phys., 128 (2003), doi 10.1063/1.2815812
 [2] Henkelmann, Jonsson, J. Chem. Phys., 111 (1999), doi 10.1063/1.480097

*/
#include <fstream>
#include <string>
#include <cmath>
#include "dimer.h"
#include "coords_io.h"
//#include "optimization_local.h"
#include "lbfgs.h"

dimermethod::kaestner_sherwood_dih::dimer_dih::dimer_dih(coords::Coordinates &cobj, coords::Representation_Main const &x, coords::Gradients_Main const &f) :
  e0(0.0), e1(0.0), C(0.0), tau(x.size()), x0(x), 
  x1(x.size()), f0(f), f1(x.size()), f_rot(x.size()), f_trans(x.size()),
  coordobject(cobj)
{ 
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::span (std::vector<coords::Representation_Main> const * rep_ptr)
{
  size_t const N(tau.size());
  for (size_t i(0U); i<N; ++i)
  {
    if (f0[i] > coords::main_gradient_type(0))
    {
      tau[i] = coords::angle_type::from_rad(coords::main_gradient_type(0) / f0[i]);
    }
  }
  //std::cout << "Span tau:" << tau << lineend;
  if (rep_ptr && (*rep_ptr).size() < tau.size())
  {
    for (auto const & tabu_direction : *rep_ptr)
    {
      tau.orthogonalize_toNormed(tabu_direction);
    }
  }
  tau.norm();
  for (size_t i(0U); i<N; ++i)
  {
    x1[i] = x0[i] + tau[i] * Config::get().dimer.distance;
    if (Config::get().general.verbosity > 99) std::cout << "x1[" << i << "] = " << x0[i] << " + " << tau[i] << "*" << Config::get().dimer.distance << " = " << x1[i] << lineend;
  }
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::update_force (coords::Representation_Main const &dihedrals, double &energy, coords::Gradients_Main &forces)
{
  coordobject.set_all_main(dihedrals);
  energy = coordobject.g();
  coordobject.to_internal();
  forces = coordobject.g_main();
}

static inline coords::float_type dot(coords::Gradients_Main const &g, coords::Representation_Main const &r)
{
  coords::Gradients_Main::size_type const N = std::min(g.size(), r.size());
  coords::float_type d = coords::float_type(0);
  for (coords::Gradients_Main::size_type i = 0; i < N; ++i)
  {
    d += g[i] * r[i].degrees();
  }
  return d;
}

static inline coords::Gradients_Main to_grad_type_from_deg(coords::Representation_Main const &r)
{
  coords::Gradients_Main::size_type const N = r.size();
  coords::Gradients_Main ret(N);
  for (coords::Gradients_Main::size_type i = 0; i < N; ++i)
  {
    ret[i] = r[i].degrees();
  }
  return ret;
}

static inline coords::Representation_Main from_grad_type_to_deg(coords::Gradients_Main const &r)
{
  coords::Representation_Main::size_type const N = r.size();
  coords::Representation_Main ret(N);
  for (coords::Gradients_Main::size_type i = 0; i < N; ++i)
  {
    ret[i] = coords::Representation_Main::value_type::from_deg(r[i]);
  }
  return ret;
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::update (bool const update_x0, bool const update_x1)
{
  if (update_x1) update_force(x1, e1, f1);
  if (update_x0) update_force(x0, e0, f0);
  tau = (x1 - x0).norm();
  size_t const N(x0.size());
  coords::Gradients_Main const f01(f1-f0), two_f01(f01*2.0);
  coords::float_type const f01_dot_tau(dot(f01, tau)), f0_dot_tau(dot(f0, tau));
  coords::Gradients_Main alt_ft2(f_trans);
  for (size_t i(0U); i<N; ++i)
  {
    f_rot[i] = -two_f01[i] + (tau[i] * 2.0*f01_dot_tau).degrees(); // ([1] Eq. (3))
    f_trans[i] = f0[i] - (tau[i] * 2.0*f0_dot_tau).degrees(); // ([1] Eq. (2))
    alt_ft2[i] -= (tau[i]*f01[i]/Config::get().dimer.distance).degrees();
  }
  C = f01_dot_tau/Config::get().dimer.distance;  // ([1] Eq. (4))
  double const factor(C > 0.0 ? (C < 1.0 ? (1.0-std::cos(PI*C))/2.0 : 1.0) : 0.0);
  //std::cout << lineend << "Update: C " << C << ", afac: " << factor << lineend;
  coords::Gradients_Main alt_ft(f_trans - to_grad_type_from_deg(tau*factor));
  if (Config::get().general.verbosity > 99U)
  {
    for (size_t i(0U); i<N; ++i)
    {
      std::cout << std::setw(18) << x0[i] << ", ";
      std::cout << std::setw(18) << x1[i] << " :: ";
      std::cout << std::setw(18) << f0[i] << ", ";
      std::cout << std::setw(18) << f1[i] << " :: ";
      std::cout << std::setw(18) << tau[i] << " :: ";
      std::cout << std::setw(18) << f_trans[i] << " :: ";
      std::cout << std::setw(18) << alt_ft[i] << " :: ";
      std::cout << std::setw(18) << alt_ft2[i] << lineend;
    }
  }
  //std::cout << std::setw(18) << "x0" << ", ";
  //std::cout << std::setw(18) << "x1" << " :: ";
  //std::cout << std::setw(18) << "f0" << ", ";
  //std::cout << std::setw(18) << "f1" << " :: ";
  //std::cout << std::setw(18) << "tau" << " :: ";
  //std::cout << std::setw(18) << "ft" << " :: ";
  //std::cout << std::setw(18) << "aft" << " :: ";
  //std::cout << std::setw(18) << "aft2" << lineend;

}

void dimermethod::kaestner_sherwood_dih::dimer_dih::rotate (coords::Representation_Main const &O, double const phi)
{
  double const cosine = std::cos(phi), sine = std::sin(phi);
  size_t const N(x0.size());
  for (size_t i(0U); i<N; ++i)
  {
    x1[i] = x0[i] + (tau[i] * cosine + O[i] * sine)*Config::get().dimer.distance; // ([2] Eq. (4))
  }
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::invert (void)
{
  size_t const M(x0.size());
  for (size_t i(0U); i<M; ++i) x1[i] = x0[i] - tau[i] * Config::get().dimer.distance;
  if (Config::get().dimer.grad_type == config::dimer::gradient_types::EXTRAPOLATE)
  {
    for (size_t i(0U); i<M; ++i) f1[i] = f0[i]+f0[i]-f1[i];
    update(false, false);
  }
  else 
  {
    update();
  }
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::move_to (coords::Representation_Main const &new_x0)
{
  size_t const M(x0.size());
  //std::cout << "Moved from " << x0 << lineend;
  for (size_t i(0U); i<M; ++i)
  {
    x0[i] = new_x0[i];
    x1[i] = x0[i] + tau[i]*Config::get().dimer.distance;
  }
  //std::cout << "Moved to " << x0 << lineend;
  update(true, true);
}

void dimermethod::kaestner_sherwood_dih::dimer_dih::advance_along_axis (coords::float_type const d)
{
  move_to(x0+(tau*d));
}

size_t dimermethod::kaestner_sherwood_dih::dimer_dih::linesearch (void)
{
  size_t i(0U);
  double step(1.0), old_frot_len(f_rot.len());
  for (; i < 50; ++i)
  {
    advance_along_axis(step);
    double const frot_len(f_rot.len());
    if (e0 > e1 || frot_len > Config::get().dimer.trans_F_rot_limit) break;
    if (frot_len < (old_frot_len*10.0) && C > 0.0) step += 1.0;
    if (Config::get().general.verbosity > 99U) std::cout << "Advance " << i << " (" << step << "): abs(f_rot): " << f_rot.len() << ", C: " << C << lineend;
    //std::cout << x0 << lineend;
    old_frot_len = frot_len;
  }
  return i;
}


dimermethod::kaestner_sherwood_dih::kaestner_sherwood_dih (coords::Coordinates &coordinates_object, 
                                                             std::vector<coords::Representation_Main> const * rep_ptr, double const) :
  m_dimer(coordinates_object, coordinates_object.main(), coordinates_object.g_main())
{
  m_dimer.span(rep_ptr);
  m_dimer.update(true);
}


bool dimermethod::kaestner_sherwood_dih::rotate_to_min_curve_cg (void)
{
  size_t const M(m_dimer.size());
  coords::Gradients_Main G_old(M), F_old(M), O_old(M);
  for (size_t i(0U); i<Config::get().dimer.maxRot; ++i)
  {
    coords::Gradients_Main G(M);
    if (i > 0U)
    {
      double const gamma = (m_dimer.f_rot-F_old).scalar(m_dimer.f_rot)/m_dimer.f_rot.scalar(m_dimer.f_rot);  // ([2] Eq. (17))
      G = m_dimer.f_rot + O_old*(gamma*G_old.len());  // ([2] Eq. (16))
    }
    else G = m_dimer.f_rot.normd();
    G.norm();
    G_old = G;
    F_old = m_dimer.f_rot;
    coords::Gradients_Main f01(m_dimer.f1 - m_dimer.f0), twof01(f01*2.0);
    double const dC_dphi(twof01.scalar(G)/Config::get().dimer.distance);  // ([1] Eq. (6))
    double const phi1(-0.5*std::atan(dC_dphi/(2.0*std::fabs(m_dimer.C)))); // ([1] Eq. (5))
    double const CONV_ANGLE(Config::get().dimer.rotationConvergence*RATIOPI180);
    if (std::fabs(phi1) < CONV_ANGLE) 
    {
      if (Config::get().general.verbosity > 9U) std::cout << "Rotation converged after " << i << " steps."  << lineend;
      return true;
    }
    else
    {
      dimer_dih dimer_tmp(m_dimer);
      dimer_tmp.rotate(from_grad_type_to_deg(G), phi1);
      dimer_tmp.update();
      double const b1(0.5*dC_dphi);  // ([1] Eq. (8))
      double const a1((m_dimer.C-dimer_tmp.C+b1*std::sin(2.0*phi1))/(1.0-std::cos(2.0*phi1)));  // ([1] Eq. (9))
      double const a0(2.0*(m_dimer.C-a1));  // ([1] Eq. (10))
      double phi(0.5*std::atan(b1/a1));
      double const C(a0/2.0 + a1*std::cos(2.0*phi) + b1*std::sin(2.0*phi)); // ([1] Eq. (7))
      if (C > m_dimer.C) phi += PI*0.5;
      if (std::fabs(phi) < CONV_ANGLE) 
      {
        if (Config::get().general.verbosity > 9U) std::cout << "Rotation converged after " << i << " steps."  << lineend;
        return true;
      }
      else
      {
        if (Config::get().general.verbosity > 49U) std::cout << "Rotating dimer (" << i << "/" << Config::get().dimer.maxRot << "; " << phi*RATIO180PI << " degrees) -> " << m_dimer.C << lineend;
        m_dimer.rotate(from_grad_type_to_deg(G), phi);
        if (Config::get().dimer.grad_type == config::dimer::gradient_types::EXTRAPOLATE || i%5 == 0)
        {
          double const sin_phi1(std::sin(phi1)), sin_phi(std::sin(phi)); // ([1] Eq. (12))
          double const f1f(std::sin(phi1-phi)/sin_phi1), ftmpf(sin_phi/sin_phi1); // ([1] Eq. (12))
          double const f0f(1.0 - std::cos(phi) - sin_phi*std::tan(phi1/2.0)); // ([1] Eq. (12))
          for (size_t j(0U); j<M; ++j)
          {
            m_dimer.f1[j] = f1f*m_dimer.f1[j] + ftmpf*dimer_tmp.f1[j] + f0f*m_dimer.f0[j];
          }
          m_dimer.update(false, false);
        }
        else
        {
          m_dimer.update();
        }
        O_old = G.orthogonalize_toNormed(to_grad_type_from_deg(m_dimer.tau));
      }
    }
  }
  return false;
}

struct Dim_Ft_Callback
{
  float E; 
  Dim_Ft_Callback (void * dimer_ptr, scon::nvect<float> const &new_x, 
              scon::nvect<float> &trans_force, size_t const, bool & go_on)
    : E(0.0)
  {
    if (dimer_ptr)
    {
      dimermethod::kaestner_sherwood_dih::dimer_dih & d_obj(*(static_cast<dimermethod::kaestner_sherwood_dih::dimer_dih*>(dimer_ptr)));
      d_obj.move_to(from_grad_type_to_deg(coords::Gradients_Main(new_x.begin(), new_x.end())));
      size_t const N(new_x.size());
      for (size_t i(0u); i<N; ++i)
      {
        trans_force[i] = static_cast<float>(d_obj.f_trans[i]);
      }
      E = -static_cast<float>(d_obj.e0);
      go_on = true;
    }
    else throw std::runtime_error("Dimer object missing for lbfgs_dimer_dih_trans_force object instantiation.");
  }
  operator float (void) const { return E; }
};

coords::Gradients_Main static inline trans_force(coords::Gradients_Main force, coords::Representation_Main tau)
{
  size_t const N(force.size());
  coords::float_type f_dot_tau(dot(force,tau));
  coords::Gradients_Main ft(N);
  for (size_t i(0U); i<N; ++i)
  {
    ft[i] = force[i] - (tau[i] * 2.0*f_dot_tau).degrees(); // ([1] Eq. (2))
  }
  return ft;
}

coords::Representation_Main dimermethod::kaestner_sherwood_dih::translate_across_transition_lbfgs (void)
{
  coords::Representation_Main const x0_start(m_dimer.x0);
  //size_t const N(m_dimer.x0.size());
  //std::cout << "Start x " << m_dimer.x0;
  if(Config::get().general.verbosity > 99U) std::cout << "Start tau: " << m_dimer.tau;
  typedef optimization::local::Wrapper<scon::nvect<float>, float, Dim_Ft_Callback> Wrapper_T;
  typedef optimization::local::linesearch::more_thuente<Wrapper_T> LineSearch_T;
  typedef optimization::local::lbfgs<LineSearch_T> Optimizer_T;
  rotate_to_min_curve_cg();
  // calculate convergence for translation from x0 and x1
  // convergence = |F|/|X|
  coords::float_type const 
    conv_x0(trans_force(m_dimer.f0, m_dimer.tau).len()/m_dimer.x0.len()), 
    conv_x1(trans_force(m_dimer.f1, m_dimer.tau).len()/m_dimer.x0.len());
  //std::cout << std::setprecision(8) << "cx0 " << conv_x0 << ", cx1 " << conv_x1 << lineend;
  // Calculate variation of convergence along dimer direction = delta_CONV/distance
  coords::float_type conv_grad((conv_x1-conv_x0)/Config::get().dimer.distance);
  //std::cout << "cgrad " << conv_grad << lineend;
  // Calculate required convergence change to leave "converged" area
  // converged if CONV < Config::get().optimization.local.bfgs.grad
  double const dconv(Config::get().optimization.local.bfgs.grad-conv_x0);
  //std::cout << "dconv " << dconv << lineend;
  // calculate required movement as 101% of distance 
  // required to raise the gradients to a non-converged value
  double const required_dx(std::sqrt(1.01*std::fabs(dconv) / std::fabs(conv_grad)));
  //std::cout << "Required movement is " << required_dx << " to leave " << Config::get().optimization.local.bfgs.grad << " convergence." << lineend;

  m_dimer.advance_along_axis(required_dx);

  //std::cout << "New TF conv " << (m_dimer.f_trans.len()/m_dimer.x0.len()) << " of " << Config::get().optimization.local.bfgs.grad << " reached." << lineend;
  
  
  Optimizer_T optimizer(&m_dimer, to_grad_type_from_deg(m_dimer.x0), m_dimer.x0.size());
  optimizer.config.epsilon = static_cast<Optimizer_T::float_t>(Config::get().optimization.local.bfgs.grad);
  optimizer.config.max_iterations = Config::get().dimer.maxStep;

  //std::cout << "xg.x: " << optimizer.xg.x;
  //std::cout << "xg.g: " << optimizer.xg.g;
  //std::cout << "xg.f: " << optimizer.xg.f << lineend;
  //std::cout << "xg.xnorm: " << optimizer.xg.xnorm << lineend;
  //std::cout << "xg.gnorm: " << optimizer.xg.gnorm << lineend;
  //std::cout << "Tau 0: " << m_dimer.tau;
  //std::cout << "Direction 0: " << optimizer.d; 
  //std::cout << "Force 0: " << m_dimer.f_trans;
  //std::cout << "f0 0: " << m_dimer.f0;
  //std::cout << "f1 0: " << m_dimer.f1;
  //std::cout << "C 0: " << m_dimer.C << lineend;
  //std::cout << "e0 0: " << m_dimer.e0 << lineend;
  //std::cout << "e1 0: " << m_dimer.e1 << lineend;

  if (m_dimer.e0 > m_dimer.e1) return coords::Representation_Main(x0_start.size());
  coords::Representation_Main x0_old(m_dimer.x0);
  
  size_t trans_iter(0U);

  for (size_t i(1U); i<=Config::get().dimer.maxStep; ++i)
  {
    //std::cout << "Tau 1a: " << m_dimer.tau;
    //std::cout << "Direction 1a: " << optimizer.d; 
    //std::cout << "Force 1a: " << m_dimer.f_trans;
    optimizer.store_xg();
    //printf("A %.3d XN %.6f GN %.6f Step %.6f\n", i, xg.xnorm, xg.gnorm, step);
    trans_iter += m_dimer.linesearch();
    //LineSearch_T::status const ls(optimizer.ls(i));
    //if (ls != LineSearch_T::status::SUCCESS)
    //{
    //  std::cout << "LS failed and is:" << ls << lineend;
    //  optimizer.xg.x = optimizer.xg_p.x;
    //  optimizer.xg.g = optimizer.xg_p.g;
    //  break;
    //}
    
    //std::cout << "Tau 1b: " << m_dimer.tau;
    //std::cout << "Direction 1b: " << optimizer.d; 
    //std::cout << "Force 1b: " << m_dimer.f_trans;
    //std::cout << "Step 1b: " << optimizer.step << lineend;
    //m_dimer.move_to(optimizer.xg.x);
    rotate_to_min_curve_cg();
    if (i > 0U && (m_dimer.x0 - x0_start).len() > (m_dimer.x1 - x0_start).len())
    {
      if (Config::get().general.verbosity > 9U) std::cout << "Inverting dimer (" << std::setprecision(10) << (m_dimer.x0 - x0_start).len() << " > " << (m_dimer.x1 - x0_start).len() << ")." << lineend;
      m_dimer.invert();
    }
    x0_old = m_dimer.x0;

    {
      coords::Gradients_Main gm(to_grad_type_from_deg(m_dimer.x0));
      optimizer.xg.update(gm, m_dimer.f_trans, static_cast<Optimizer_T::float_t>(-m_dimer.e0));
    }
    
    if (i > 0U && m_dimer.e0 > m_dimer.e1)
    {
      //m_dimer.advance_along_axis(1.0);
      m_dimer.coordobject.set_all_main(m_dimer.x1);
      m_dimer.coordobject.g();
      break;
    }
    //if (optimizer.convergence()) break;
    optimizer.assignHessianCorrection();
    optimizer.updateHessian(i);
    //std::cout << "Tau 2: " << m_dimer.tau;
    //std::cout << "Direction 2: " << optimizer.d; 
    //std::cout << "Force 2: " << m_dimer.f_trans;
    optimizer.step = 1.0;
  }
  return (m_dimer.x0-x0_start);
}

void dimermethod::kaestner_sherwood_dih::translate_across_transition_cg (void)
{
  size_t const M(m_dimer.size());
  coords::Gradients_Main  G_old(M), F_old(M);
  coords::Representation_Main x0_old(m_dimer.x0);
  for (size_t i(0U); i<Config::get().dimer.maxStep; ++i)
  {
    rotate_to_min_curve_cg();
    if (i > 0U && (m_dimer.x0-x0_old).len() > (m_dimer.x1-x0_old).len())
    {
      if (Config::get().general.verbosity > 9U) std::cout << "Inverting dimer." << lineend;
      m_dimer.invert();
    }
    double const gamma(i%10==0U ? 0.0 : std::max(0.0, ((m_dimer.f_trans-F_old).scalar(m_dimer.f_trans)/m_dimer.f_trans.scalar(m_dimer.f_trans))));
    coords::Gradients_Main G((m_dimer.f_trans + G_old*gamma));
    F_old = m_dimer.f_trans;
    G_old = G;
    G.set_len(Config::get().dimer.distance);
  }
}
