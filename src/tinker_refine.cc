#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include "Scon/scon_linkedcell.h"
#include "tinker_refine.h"
#include "tinker_parameters.h"
#include "coords.h"
#include "configuration.h"
#include "Scon/scon_chrono.h"

//void tinker::refine::refined::refine(coords::Coordinates const &cobj,
//  tinker::parameter::parameters const &pobj)
//{
//  coords = &cobj;
//  params = &pobj;
//  // clear refined data
//  clear();
//  // get nonbonded parameter matrices
//  m_vdwc_matrices = pobj.vdwc_matrices();
//  //std::cout<<"REFINE"<<'\n';
//  std::size_t const N_a(cobj.atoms().size());
//  m_removes.resize(N_a);
//  m_red_types.resize(N_a);
//  for (std::size_t i(0u); i < 5u; ++i)
//  {
//    m_relations[i].resize(N_a);
//  }
//
//  /* ------------
//    a
//    -> b
//       -> c
//          -> d
//             -> e
//  ------------- */
//
//  for (std::size_t a(0u); a < N_a; ++a)
//  {
//    // stuff with a
//    m_red_types[a] = pobj.contract_type(cobj.atoms(a).energy_type());
//    std::size_t N_b(cobj.atoms(a).bonds().size());
//    if (N_b == 3) refine_imptor_improper_of_atom(a);
//    for (std::size_t b1(0u); b1 < N_b; ++b1)
//    {
//      std::size_t const b(cobj.atoms(a).bonds(b1));
//      // stuff with b
//      if (b > a)
//      {
//        add_relation(a, b, R12);
//        if (!find_bond(a, b))
//        {
//          std::cout << "No bond parameters found for [Number:" << a + 1 << "," << b + 1 << "]";
//          std::cout << "[Type:" << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type() << "]";
//          std::cout << "[Type/Class:" << params->type(cobj.atoms(a).energy_type(), tinker::BOND, true);
//          std::cout << "," << params->type(cobj.atoms(b).energy_type(), tinker::BOND, true) << "]\n";
//        }
//      }
//      std::size_t N_c(cobj.atoms(b).bonds().size());
//      for (std::size_t b2(0u); b2 < N_c; ++b2)
//      {
//        std::size_t const c(cobj.atoms(b).bonds(b2));
//        // stuff with c
//        if (a == c) continue;
//        if (c > a)
//        {
//          add_relation(a, c, R13);
//          if (!cobj.atoms().sub_io() || !cobj.atoms().sub_io_transition(a, c))
//          {
//            if (!find_angle(a, b, c))
//            {
//              std::cout << "No angle parameters found for [Number:" << a + 1 << "," << b + 1 << "," << c + 1 << "]";
//              std::cout << "[Type:" << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type() << "," << cobj.atoms(c).energy_type() << "]";
//              std::cout << "[Type/Class:" << params->type(cobj.atoms(a).energy_type(), tinker::ANGLE, true);
//              std::cout << "," << params->type(cobj.atoms(b).energy_type(), tinker::ANGLE, true);
//              std::cout << "," << params->type(cobj.atoms(c).energy_type(), tinker::ANGLE, true) << "]\n";
//            }
//            find_urey(a, b, c);
//            find_strbend(a, b, c);
//          }
//
//        }
//        std::size_t N_d(cobj.atoms(c).bonds().size());
//        for (std::size_t b3(0u); b3 < N_d; ++b3)
//        {
//          std::size_t const d(cobj.atoms(c).bonds(b3));
//          // stuff with d
//          if (d == b || d == a) continue;
//          if (d > a)
//          {
//            add_relation(a, d, R14);
//            if (!cobj.atoms().sub_io() || !cobj.atoms().sub_io_transition(a, d))
//            {
//              if (!find_torsion(a, b, c, d))
//              {
//                std::cout << "No torsion parameters found for [Number:" << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << "]";
//                std::cout << "[Types: " << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type();
//                std::cout << "," << cobj.atoms(c).energy_type() << "," << cobj.atoms(d).energy_type() << "]";
//                std::cout << "[Type/Class: " << params->type(cobj.atoms(a).energy_type(), tinker::TORSION, true);
//                std::cout << "," << params->type(cobj.atoms(b).energy_type(), tinker::TORSION, true);
//                std::cout << "," << params->type(cobj.atoms(c).energy_type(), tinker::TORSION, true);
//                std::cout << "," << params->type(cobj.atoms(d).energy_type(), tinker::TORSION, true) << "]\n";
//              }
//            }
//            // todo find_opbend()
//            if (!find_opbend(a, b, c, d))
//            {
//              //std::cout << "No oopb parameters found for [Number:" << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << "]"<<'\n';
//
//            }
//          }
//          std::size_t N_e(cobj.atoms(d).bonds().size());
//          for (std::size_t b4(0u); b4 < N_e; ++b4)
//          {
//            std::size_t const e(cobj.atoms(d).bonds(b4));
//            // stuff with e
//            if (e == c || e == b || e == a) continue;
//            if (e > a) add_relation(a, e, R15);
//          }
//        }
//      }
//    }
//  }
//
//  if (!params->multipoles().empty()) refine_mp();
//  if (!params->polarizes().empty()) refine_pol();
//  refine_nb(cobj);
//}

tinker::refine::refined::refined(coords::Coordinates const & cobj, tinker::parameter::parameters const & pobj)
  : m_cparams(pobj)
{

  // get nonbonded parameter matrices
  this->m_vdwc_matrices = pobj.vdwc_matrices();
  //std::cout<<"REFINE"<<'\n';
  std::size_t const N_atoms(cobj.atoms().size());

  using namespace tinker::refine;

  // Refined relations 
  // 5u -> (-11- -12- -13- -14- -15-)
  for (std::size_t i(0u); i < 5u; ++i)
  {
    this->m_relations[i].resize(N_atoms);
  }

  // removed relations
  //vector_size_2d m_removes;
  this->m_removes.resize(N_atoms);
  //vector_size_1d m_red_types;
  this->m_red_types.resize(N_atoms);


  /* ------------
    a
    -> b
       -> c
          -> d
             -> e
  ------------- */

  for (std::size_t a(0u); a < N_atoms; ++a)
  {
    // stuff with a
    this->m_red_types[a] = pobj.contract_type(cobj.atoms(a).energy_type());
    std::size_t N_b(cobj.atoms(a).bonds().size());
    if (N_b == 3)
    {
      std::pair<vector_imptors, vector_improper> returns = refine_imptor_improper_of_atom(cobj, pobj, a);
      this->m_impropers.insert(std::end(m_impropers), std::begin(std::get< vector_improper>(returns)), std::end(std::get< vector_improper>(returns)));
      this->m_imptors.insert(std::end(m_imptors), std::begin(std::get< vector_imptors>(returns)), std::end(std::get< vector_imptors>(returns)));
    }

    for (std::size_t b1(0u); b1 < N_b; ++b1)
    {
      std::size_t const b(cobj.atoms(a).bonds(b1));
      // stuff with b
      if (b > a)
      {
        this->add_relation(pobj, a, b, R12);
        //if (!find_bond(a, b))
        auto bondsFound = find_bond(cobj, pobj, a, b);
        if(bondsFound == boost::none)
        {
          std::cout << "No bond parameters found for [Number:" << a + 1 << "," << b + 1 << "]";
          std::cout << "[Type:" << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type() << "]";
          std::cout << "[Type/Class:" << pobj.type(cobj.atoms(a).energy_type(), tinker::BOND, true);
          std::cout << "," << pobj.type(cobj.atoms(b).energy_type(), tinker::BOND, true) << "]\n";
        }
        else
        {
          this->m_bonds.emplace_back(bondsFound.get());
        }
      }
      std::size_t N_c(cobj.atoms(b).bonds().size());
      for (std::size_t b2(0u); b2 < N_c; ++b2)
      {
        std::size_t const c(cobj.atoms(b).bonds(b2));
        // stuff with c
        if (a == c) continue;
        if (c > a)
        {
          this->add_relation(pobj, a, c, R13);
          if (!cobj.atoms().sub_io() || !cobj.atoms().sub_io_transition(a, c))
          {
            auto anglesFound = find_angle(cobj, pobj, a, b,c);
            if (anglesFound == boost::none)
            {
              std::cout << "No angle parameters found for [Number:" << a + 1 << "," << b + 1 << "," << c + 1 << "]";
              std::cout << "[Type:" << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type() << "," << cobj.atoms(c).energy_type() << "]";
              std::cout << "[Type/Class:" << pobj.type(cobj.atoms(a).energy_type(), tinker::ANGLE, true);
              std::cout << "," << pobj.type(cobj.atoms(b).energy_type(), tinker::ANGLE, true);
              std::cout << "," << pobj.type(cobj.atoms(c).energy_type(), tinker::ANGLE, true) << "]\n";
            }
            else
            {
              this->m_angles.emplace_back(anglesFound.get());
            }
            auto ureyFound = find_urey(cobj, pobj, a, b, c);
            if (ureyFound != boost::none)
            {
              m_ureys.emplace_back(ureyFound.get());
            }
            auto strbendFound = find_strbend(cobj, pobj, a, b, c);
            if (strbendFound != boost::none)
            {
              m_strbends.emplace_back(strbendFound.get());
            }
          }

        }
        std::size_t N_d(cobj.atoms(c).bonds().size());
        for (std::size_t b3(0u); b3 < N_d; ++b3)
        {
          std::size_t const d(cobj.atoms(c).bonds(b3));
          // stuff with d
          if (d == b || d == a) continue;
          if (d > a)
          {
            this->add_relation(pobj, a, d, R14);
            if (!cobj.atoms().sub_io() || !cobj.atoms().sub_io_transition(a, d))
            {
              auto torsionFound = find_torsion(cobj, pobj, a, b, c,d);
              if (torsionFound == boost::none)
              {
                std::cout << "No torsion parameters (or parameters exactly equal to zero) found for [Number:" << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << "]";
                std::cout << "[Types: " << cobj.atoms(a).energy_type() << "," << cobj.atoms(b).energy_type();
                std::cout << "," << cobj.atoms(c).energy_type() << "," << cobj.atoms(d).energy_type() << "]";
                std::cout << "[Type/Class: " << pobj.type(cobj.atoms(a).energy_type(), tinker::TORSION, true);
                std::cout << "," << pobj.type(cobj.atoms(b).energy_type(), tinker::TORSION, true);
                std::cout << "," << pobj.type(cobj.atoms(c).energy_type(), tinker::TORSION, true);
                std::cout << "," << pobj.type(cobj.atoms(d).energy_type(), tinker::TORSION, true) << "]\n";
              }
              else
              {
                m_torsions.emplace_back(torsionFound.get());
              }
            }
            // todo find_opbend()
            auto opbendFound = find_opbend(cobj, pobj, a, b, c, d);
            if (opbendFound != boost::none)
            {
              this->m_opbends.emplace_back(opbendFound.get());
            }
          }
          std::size_t N_e(cobj.atoms(d).bonds().size());
          for (std::size_t b4(0u); b4 < N_e; ++b4)
          {
            std::size_t const e(cobj.atoms(d).bonds(b4));
            // stuff with e
            if (e == c || e == b || e == a) continue;
            if (e > a) this->add_relation(pobj, a, e, R15);
          }
        }
      }
    }
  }

  if (!pobj.multipoles().empty()) 
    this->m_multipole_vec.emplace_back(refine_mp(cobj,pobj));
  if (!pobj.polarizes().empty()) 
    this->m_polarize_vec.emplace_back(refine_pol(cobj, pobj));
  refine_nb(cobj);
}

void tinker::refine::refined::add_relation(tinker::parameter::parameters const & pobj, std::size_t const atom,
  std::size_t const related, std::size_t const relation)
{
  if (!has_tighter_relation(atom, related, relation))
  {  // check whether we have a tighter relation than the specified one betweem atom and related
    scon::sorted::insert_unique(m_relations[relation][atom], related);
    if (!pobj.vdwc_used(relation))
    {
      scon::sorted::insert_unique(m_removes[atom], related);
    }
    // make sure this relation is the tightest
    remove_loose_relations(atom, related, relation);
  }
}

bool tinker::refine::refined::has_tighter_relation(std::size_t const atom,
  std::size_t const related, std::size_t const relation_to_check)
{
  for (std::size_t i(0u); i < relation_to_check; ++i)
  {
    if (scon::sorted::exists(m_relations[i][atom], related)) return true;
  }
  return false;
}

void tinker::refine::refined::remove_loose_relations(std::size_t const atom,
  std::size_t const related, std::size_t const relation_to_check)
{
  for (std::size_t i(relation_to_check + 1); i < m_relations.size(); ++i)
  {
    if (!m_relations[i][atom].empty())
    {
      if ((m_relations[i][atom].back() < related) || (related < m_relations[i][atom].front())) continue;
      std::size_t find_index = scon::sorted::find(m_relations[i][atom], related);
      if (find_index == m_relations[i][atom].size() || related < m_relations[i][atom][find_index]) continue;
      m_relations[i][atom].erase(m_relations[i][atom].begin() + static_cast<std::ptrdiff_t>(find_index));
    }
  }
}

//ok
tinker::refine::vector_improper tinker::refine::find_improper(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t center, std::size_t a, std::size_t b, std::size_t c)
{
  std::size_t const tc(params.type(coords.atoms(center).energy_type(), tinker::IMPROPER)),
    t1(params.type(coords.atoms(a).energy_type(), tinker::IMPROPER)),
    t2(params.type(coords.atoms(b).energy_type(), tinker::IMPROPER)),
    t3(params.type(coords.atoms(c).energy_type(), tinker::IMPROPER));
  std::size_t sym(0u);
  tinker::refine::vector_improper m_impropers;
  using namespace tinker::refine;

  for (auto const & p : params.impropers())
  {
    if (p.center == tc)
    {
      double dsym(1.0);
      if (t1 == t2 || t1 == t3 || t2 == t3) dsym = 2.0;
      if (t1 == t2 && t1 == t3 && t2 == t3) dsym = 6.0;
      if (p.check_lig(t1, t2, t3))
      {
        m_impropers.push_back(types::improper(center, a, b, c, p, dsym));
        ++sym;
      }
      if (p.check_lig(t1, t3, t2))
      {
        m_impropers.push_back(types::improper(center, a, c, b, p, dsym));
        ++sym;
      }
      if (p.check_lig(t2, t1, t3))
      {
        m_impropers.push_back(types::improper(center, b, a, c, p, dsym));
        ++sym;
      }
      if (p.check_lig(t2, t3, t1))
      {
        m_impropers.push_back(types::improper(center, b, c, a, p, dsym));
        ++sym;
      }
      if (p.check_lig(t3, t2, t1))
      {
        m_impropers.push_back(types::improper(center, c, b, a, p, dsym));
        ++sym;
      }
      if (p.check_lig(t3, t1, t2))
      {
        m_impropers.push_back(types::improper(center, c, a, b, p, dsym));
        ++sym;
      }
    }
  }

  if (sym < 1)
  {
    for (auto const & p : params.impropers())
    {
      if (p.center == tc && p.ligand[0] == 0 && p.ligand[1] == 0)
      {
        if (p.ligand[2] == t1 || p.ligand[2] == t2 || p.ligand[2] == t3)
        {
          m_impropers.push_back(types::improper(center, a, b, c, p, 3.0));
          m_impropers.push_back(types::improper(center, b, c, a, p, 3.0));
          m_impropers.push_back(types::improper(center, c, a, b, p, 3.0));
          sym += 3;
        }
      }
    }
  }

  if (sym < 1)
  {
    for (auto const & p : params.impropers())
    {
      if (p.center == tc && p.ligand[0] == 0 && p.ligand[1] == 0 && p.ligand[2] == 0)
      {
        m_impropers.push_back(types::improper(center, a, b, c, p, 3.0));
        m_impropers.push_back(types::improper(center, b, c, a, p, 3.0));
        m_impropers.push_back(types::improper(center, c, a, b, p, 3.0));
        sym += 3;
      }
    }
  }

  return m_impropers;
}
//ok
boost::optional<tinker::refine::types::binary_quadratic> tinker::refine::find_bond(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b)
  {
    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::BOND)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::BOND));

    tinker::refine::types::binary_quadratic pot;
    pot.atoms[0] = a;
    pot.atoms[1] = b;
    for (auto const & i : params.bonds())
    {
      if (i.check(ta, tb))
      {
        pot.force = i.f;
        pot.ideal = i.ideal;
        return pot;
      }
    }
    return boost::none;
  }
//ok
boost::optional<tinker::refine::types::ternary_quadratic> tinker::refine::find_angle(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c)
  {
    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::ANGLE)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::ANGLE)),
      tc(params.type(coords.atoms(c).energy_type(), tinker::ANGLE));
    tinker::refine::types::ternary_quadratic pot;
    pot.atoms[0] = a;
    pot.atoms[1] = b;
    pot.atoms[2] = c;
    for (auto const & i : params.angles())
    {
      if (i.check(ta, tb, tc))
      {
        pot.force = i.f;
        pot.ideal = i.ideal;
        return pot;
      }
    }
    return boost::none;
  }
//ok
tinker::refine::vector_imptors tinker::refine::find_imptors(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t center, std::size_t a, std::size_t b, std::size_t c)
{
  std::size_t const tc(params.type(coords.atoms(center).energy_type(), tinker::IMPTORS)),
    t1(params.type(coords.atoms(a).energy_type(), tinker::IMPTORS)),
    t2(params.type(coords.atoms(b).energy_type(), tinker::IMPTORS)),
    t3(params.type(coords.atoms(c).energy_type(), tinker::IMPTORS));
  std::size_t sym(0u);
  using namespace tinker::refine;
  tinker::refine::vector_imptors m_imptors;

  for (auto const & p : params.imptors())
  {
    if (p.center == tc)
    {
      double dsym(1.0);
      if (t1 == t2 || t1 == t3 || t2 == t3) dsym = 2.0;
      if (t1 == t2 && t1 == t3 && t2 == t3) dsym = 6.0;
      if (p.check_lig(t1, t2, t3))
      {
        m_imptors.push_back(types::imptor(center, a, b, c, p, dsym));
        ++sym;
      }
      if (p.check_lig(t1, t3, t2))
      {
        m_imptors.push_back(types::imptor(center, a, c, b, p, dsym));
        ++sym;
      }
      if (p.check_lig(t2, t1, t3))
      {
        m_imptors.push_back(types::imptor(center, b, a, c, p, dsym));
        ++sym;
      }
      if (p.check_lig(t2, t3, t1))
      {
        m_imptors.push_back(types::imptor(center, b, c, a, p, dsym));
        ++sym;
      }
      if (p.check_lig(t3, t2, t1))
      {
        m_imptors.push_back(types::imptor(center, c, b, a, p, dsym));
        ++sym;
      }
      if (p.check_lig(t3, t1, t2))
      {
        m_imptors.push_back(types::imptor(center, c, a, b, p, dsym));
        ++sym;
      }
    }
  }

  if (sym < 1)
  {
    for (auto const & p : params.imptors())
    {
      if (p.center == tc && p.ligand[0] == 0 && p.ligand[1] == 0)
      {
        if (p.ligand[2] == t1 || p.ligand[2] == t2 || p.ligand[2] == t3)
        {
          m_imptors.push_back(types::imptor(center, a, b, c, p, 3.0));
          m_imptors.push_back(types::imptor(center, b, c, a, p, 3.0));
          m_imptors.push_back(types::imptor(center, c, a, b, p, 3.0));
          sym += 3;
        }
      }
    }
  }

  if (sym < 1)
  {
    for (auto const & p : params.imptors())
    {
      if (p.center == tc && p.ligand[0] == 0 && p.ligand[1] == 0 && p.ligand[2] == 0)
      {
        m_imptors.push_back(types::imptor(center, a, b, c, p, 3.0));
        m_imptors.push_back(types::imptor(center, b, c, a, p, 3.0));
        m_imptors.push_back(types::imptor(center, c, a, b, p, 3.0));
        sym += 3;
      }
    }
  }

  return m_imptors;
}
//ok
boost::optional<tinker::refine::types::torsion> tinker::refine::find_torsion(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c, std::size_t d)
  {
    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::TORSION)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::TORSION)),
      tc(params.type(coords.atoms(c).energy_type(), tinker::TORSION)),
      td(params.type(coords.atoms(d).energy_type(), tinker::TORSION));

    using namespace tinker::refine;
    types::torsion pot;
    pot.atoms[0] = a;
    pot.atoms[1] = b;
    pot.atoms[2] = c;
    pot.atoms[3] = d;
    std::size_t maxcheck(0u);
    for (auto const & i : params.torsions())
    {
      std::size_t const check(i.check(ta, tb, tc, td));
      if (check > maxcheck)
      {
        pot.p = i;
        maxcheck = check;
        if (check == 4)  // only "real" atom type parameters
        {
          if (!i.empty()) return pot;
        }
      }
    }
    if (maxcheck > 0U) // not only "real" atom type parameter but least number of 0s in torsion definition
    {
      if (!pot.p.empty()) return pot;
    }
    return boost::none;
  }
//ok
boost::optional<tinker::refine::types::opbend> tinker::refine::find_opbend(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c, std::size_t d)
  {
    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::OPBEND)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::OPBEND))
      //, tc(params.type(coords.atoms(c).energy_type(), tinker::OPBEND))
      //, td(params.type(coords.atoms(d).energy_type(), tinker::OPBEND))
      ;
    using namespace tinker::refine;

    types::opbend pot;
    pot.atoms[0] = a;
    pot.atoms[1] = b;
    pot.atoms[2] = c;
    pot.atoms[3] = d;

    for (auto const & i : params.opbends())
    {
      if (i.check(ta, tb))
      {
        pot.force = i.f;
        return pot;
      }
      else
      {
        return boost::none;
      }
    }

    return boost::none;
  }
//ok
boost::optional<tinker::refine::types::strbend> tinker::refine::find_strbend(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c)
  {
    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::STRBEND)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::STRBEND)),
      tc(params.type(coords.atoms(c).energy_type(), tinker::STRBEND));
    using namespace tinker::refine;

    types::strbend pot;
    pot.atoms[0] = a;
    pot.atoms[1] = b;
    pot.atoms[2] = c;
    for (auto const & i : params.strbends())
    {
      if (i.check(ta, tb, tc) == 0U)
      {
        pot.force = i.f;
        pot.ff = i.ff;
        return pot;
      }
      else if (i.check(ta, tb, tc) == 1U)
      {
        pot.force = i.ff;
        pot.ff = i.f;
        return pot;
      }
    }
    return boost::none;
  }
//ok
boost::optional<tinker::refine::types::binary_quadratic> tinker::refine::find_urey(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t a, std::size_t b, std::size_t c)
  {
    // check if a-->c (c-->a) is in->out or out->in transition
    if (coords.atoms().sub_io() && coords.atoms().sub_io_transition(a, c))
    {
      return boost::none;
    }
    using namespace tinker::refine;

    std::size_t const ta(params.type(coords.atoms(a).energy_type(), tinker::UREYBRAD)),
      tb(params.type(coords.atoms(b).energy_type(), tinker::UREYBRAD)),
      tc(params.type(coords.atoms(c).energy_type(), tinker::UREYBRAD));
    types::binary_quadratic pot;
    pot.atoms[0] = a;
    pot.atoms[1] = c;
    for (auto const & i : params.ureybrads())
    {
      if (i.check(ta, tb, tc))
      {
        pot.force = i.f;
        pot.ideal = i.ideal;
        return pot;
      }
    }
    return boost::none;
  }

std::pair<tinker::refine::vector_imptors, tinker::refine::vector_improper> tinker::refine::refine_imptor_improper_of_atom(coords::Coordinates const & coords, tinker::parameter::parameters const & params, std::size_t atom)
{
  std::size_t const a(atom), b(coords.atoms(atom).bonds(0)),
    c(coords.atoms(atom).bonds(1)), d(coords.atoms(atom).bonds(2));
  tinker::refine::vector_improper m_impropers = find_improper(coords, params, a, b, c, d);
  tinker::refine::vector_imptors m_imptors = find_imptors(coords, params, a, b, c, d);
  return std::make_pair(m_imptors, m_impropers);
}

void tinker::refine::refined::refine_nb(coords::Coordinates const & coords)
{
  //scon::chrono::high_resolution_timer rnbt;

  if (m_cparams.vdwc_used(R12)) build_pairs_direct<R12>(coords);
  else if (m_cparams.vdwc_used(R13)) build_pairs_direct<R13>(coords);
  else if (m_cparams.vdwc_used(R14)) build_pairs_direct<R14>(coords);
  else if (m_cparams.vdwc_used(R15)) build_pairs_direct<R15>(coords);
  else build_pairs_direct<R1N>(coords);


  //std::cout << "refine_nb time: " << rnbt << "\n";
  // Todo: build pairs using new-linkedcells (todo: new linkedcells)
}

template<tinker::refine::refined::rel RELATION>
bool tinker::refine::refined::add_pair(coords::Coordinates const & coords, std::size_t const row, std::size_t const col, std::array<std::size_t, 5u> const & to_matrix_id)
{
  ::tinker::parameter::parameters const* params = &m_cparams;
  if (
      (!coords.atoms().sub_io() || !coords.atoms().sub_io_transition(row, col))
      && !scon::sorted::exists(m_removes[col], row)
      && !(
        Config::get().energy.remove_fixed
        && coords.atoms(row).fixed()
        && coords.atoms(col).fixed()
        )
      )
  {
    types::nbpair pair(row, col);
    if (RELATION == R12 && to_matrix_id[R12] > 0 && scon::sorted::exists(m_relations[R12][col], row))
    {
      if (params->vdwc_used(R12))
      {
        m_pair_matrices[to_matrix_id[R12]].pair_matrix(coords.atoms(row).system(), coords.atoms(col).system()).push_back(pair);
      }
    }
    else if ((RELATION == R12 || RELATION == R13) && to_matrix_id[R13] > 0 && scon::sorted::exists(m_relations[R13][col], row))
    {
      if (params->vdwc_used(R13))
      {
        m_pair_matrices[to_matrix_id[R13]].pair_matrix(coords.atoms(row).system(), coords.atoms(col).system()).push_back(pair);
      }
    }
    else if ((RELATION == R12 || RELATION == R13 || RELATION == R14) && to_matrix_id[R14] > 0 && scon::sorted::exists(m_relations[R14][col], row))
    {
      if (params->vdwc_used(R14))
      {
        m_pair_matrices[to_matrix_id[R14]].pair_matrix(coords.atoms(row).system(), coords.atoms(col).system()).push_back(pair);
      }
    }
    else if ((RELATION == R12 || RELATION == R13 || RELATION == R14 || RELATION == R15) && to_matrix_id[R15] > 0 && scon::sorted::exists(m_relations[R15][col], row))
    {
      if (params->vdwc_used(R15))
      {
        m_pair_matrices[to_matrix_id[R15]].pair_matrix(coords.atoms(row).system(), coords.atoms(col).system()).push_back(pair);
      }
    }
    else
    {
      m_pair_matrices[0u].pair_matrix(coords.atoms(row).system(), coords.atoms(col).system()).push_back(pair);
    }
    return true;
  }
  return false;
}

template<tinker::refine::refined::rel RELATION>
void tinker::refine::refined::build_pairs_direct(coords::Coordinates const & coords)
{
  m_pair_matrices.clear();
  std::size_t const num_sys = coords.subsystems().size();
  std::size_t const num_ia = (num_sys*num_sys + num_sys) / 2;
  m_pair_matrices.push_back(types::nbpm(num_ia, 5u));
  std::array<std::size_t, 5u> to_matrix_id = { { 0u, 0u, 0u, 0u, 0u } };
  // determine the pair matrices to be applied
  for (std::size_t i(RELATION); i < 5u; ++i)
  {
    bool const bool1 = m_cparams.vdwc_used(i);
    bool const bool2 = !m_vdwc_matrices[i].empty();
    if (bool1 && bool2)
    {
      if (!m_cparams.vdwc_scaled(i) && (i != 3U || !m_cparams.has_vdw14s()))
      {
        to_matrix_id[i] = 0u;
      }
      else
      {
        bool create_matrix(true);
        std::size_t const mpms(m_pair_matrices.size());
        for (std::size_t m_id(0u); m_id < mpms; ++m_id)
        {
          if (m_pair_matrices[m_id].param_matrix_id < 5u && m_cparams.vdwc_same(i, m_pair_matrices[m_id].param_matrix_id))
          { // if the same scaling is applied as for another existing matrix we do not create a new one
            create_matrix = false;
            to_matrix_id[i] = m_id;
            break;
          }
        }
        if (create_matrix)
        {
          m_pair_matrices.push_back(types::nbpm(num_ia, i));
          to_matrix_id[i] = m_pair_matrices.size() - 1;
        }
      }
    }
  }
  for (std::size_t k = 0; k < num_sys; ++k)
  {
    for (std::size_t l = 0; l<num_sys; ++l)
    {
      if (k != l)
      {
        m_pair_matrices[0].pair_matrix(k, l).reserve(coords.subsystems(k).size()*coords.subsystems(l).size());
      }
      else
      {
        std::size_t const ks = coords.subsystems(k).size();
        m_pair_matrices[0].pair_matrix(k, l).reserve((ks*(ks - 1u)) / 2u);
      }
    }
  }

  if (Config::get().energy.cutoff > 500.0 || Config::get().periodics.periodic)   // no linked cell algorithm
  {                                                    // if cutoff > 500 or periodic boundaries activated
    const std::size_t N = coords.size(), M = (N*N - N) / 2;
    for (std::size_t i(0u), atom_row(1u), relation_col(0u); i < M; ++i)
    {
      add_pair<RELATION>(coords, atom_row, relation_col, to_matrix_id);
      ++relation_col;
      if (relation_col == atom_row)
      {
        relation_col = 0;
        ++atom_row;
      }
    }
  }
  else    // linked cell algorithm
  {
    using cells_type = scon::linked::Cells < coords::float_type, coords::Cartesian_Point, coords::Representation_3D >;
    cells_type atmcells(
      coords.xyz(),
      Config::get().energy.cutoff,
      Config::get().periodics.periodic, Config::get().periodics.pb_box, coords::float_type(0),
      scon::linked::fragmentation::half);
    std::size_t const N = coords.size();
    coords::Cartesian_Point const halfbox(Config::get().periodics.pb_box / 2.);
    //#pragma omp parallel     // to be reactivated at some point in the future?
    for (std::size_t i = 0; i < N; ++i)
    {
      auto box_of_i = atmcells.box_of_element(i);
      for (auto j : box_of_i.adjacencies())
      {
        if (j >= 0)
        {
          std::size_t const uj = static_cast<std::size_t>(j);
          if (uj < i)
          {
            /*if (Config::get().energy.periodic)
            {
              coords::Cartesian_Point d(scon::abs(coords.xyz(i) - coords.xyz(uj)));
              if (d.x() > halfbox.x()) d.x() = std::abs(d.x() - Config::get().energy.pb_box.x());
              if (d.y() > halfbox.y()) d.y() = std::abs(d.y() - Config::get().energy.pb_box.y());
              if (d.z() > halfbox.z()) d.z() = std::abs(d.z() - Config::get().energy.pb_box.z());
              if (geometric_length(d) < Config::get().energy.cutoff + 5.0)
              {
                add_pair<RELATION>(i, uj, to_matrix_id);
              }
            }
            else
            {*/
              add_pair<RELATION>(coords, i, uj, to_matrix_id);
            /*}*/
          }
        }
      }
    }
  }
}

tinker::refine::vector_multipole tinker::refine::refine_mp(coords::Coordinates const & coords, tinker::parameter::parameters const & params)
{
  tinker::refine::vector_multipole m_multipole;
  double const B = 0.52917721092000003;
  std::size_t const N(coords.atoms().size());

  refine::types::multipole refined_multipole;

  refined_multipole.npole = 0;
  std::cout << "NUMBER_ATOMS " << N << '\n';
  for (size_t a(0U); a < N; ++a)
  {
    std::size_t const type_of_a(params.type(coords.atoms(a).energy_type(), tinker::MULTIPOLE));
    refined_multipole.center = a;
    for (auto const & mpp : params.multipoles())
    {
      if (type_of_a == mpp.index[0])
      {
        if (mpp.index[1] == 0)
        {
          //refined_multipole.p_nonrot = &mpp;
          refined_multipole.axes.z() = a;
          refined_multipole.p_rot.charge = mpp.charge;
          refined_multipole.p_rot.dipole = mpp.dipole;
          refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
          refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
          refined_multipole.p_nonrot = &refined_multipole.p_rot;
          refined_multipole.p_rot.axt = mpp.axt;
          refined_multipole.npole++;
        }
        else
        {
          bool found_par_b(false);
          for (auto const b : coords.atoms(a).bonds())
          {

            std::size_t const type_of_b(params.type(coords.atoms(b).energy_type(), tinker::MULTIPOLE));
            if (type_of_b == mpp.index[1])
            {
              refined_multipole.axes.z() = b;
              if (mpp.index[2] == 0)
              {
                //refined_multipole.p_nonrot = &mpp;
                found_par_b = true;
                refined_multipole.p_rot.charge = mpp.charge;
                refined_multipole.p_rot.dipole = mpp.dipole;
                refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
                refined_multipole.p_nonrot = &refined_multipole.p_rot;
                refined_multipole.p_rot.axt = mpp.axt;
                refined_multipole.npole++;
                break;
              }
              else
              {
                // check 1-2 multipole axes 
                bool found_par_c(false);
                for (auto const c : coords.atoms(a).bonds())
                {
                  if (c != b)
                  {
                    std::size_t const type_of_c(params.type(coords.atoms(c).energy_type(), tinker::MULTIPOLE));
                    if (type_of_c == mpp.index[2])
                    {
                      refined_multipole.axes.x() = c;
                      if (mpp.index[3] == 0)
                      {
                        //refined_multipole.p_nonrot = &mpp;
                        refined_multipole.p_rot.charge = mpp.charge;
                        refined_multipole.p_rot.dipole = mpp.dipole;
                        refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                        refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
                        refined_multipole.p_nonrot = &refined_multipole.p_rot;
                        refined_multipole.p_rot.axt = mpp.axt;
                        found_par_b = found_par_c = true;
                        refined_multipole.npole++;
                        break;
                      }
                      else
                      {
                        for (auto const d : coords.atoms(a).bonds())
                        {
                          if (d != c && d != b)
                          {
                            std::size_t const type_of_d(params.type(coords.atoms(d).energy_type(), tinker::MULTIPOLE));
                            if (type_of_d == mpp.index[3])
                            {
                              refined_multipole.axes.y() = d;
                              //refined_multipole.p_nonrot = &mpp;
                              refined_multipole.p_rot.charge = mpp.charge;
                              refined_multipole.p_rot.dipole = mpp.dipole;
                              refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                              refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
                              refined_multipole.p_nonrot = &refined_multipole.p_rot;
                              refined_multipole.p_rot.axt = mpp.axt;



                              found_par_b = found_par_c = true;
                              refined_multipole.npole++;
                              break;
                            }
                          }
                        }
                        if (found_par_c) break;
                      }
                    }
                  }
                }
                // check 1-3 multipole axes if no parameter has been found
                if (!found_par_c)
                {
                  for (auto const c : coords.atoms(b).bonds())
                  {
                    if (c != a)
                    {
                      std::size_t const type_of_c(params.type(coords.atoms(c).energy_type(), tinker::MULTIPOLE));
                      if (type_of_c == mpp.index[2])
                      {
                        refined_multipole.axes.x() = c;
                        if (mpp.index[3] == 0)
                        {
                          //refined_multipole.p_nonrot = &mpp;
                          refined_multipole.p_rot.charge = mpp.charge;
                          refined_multipole.p_rot.dipole = mpp.dipole;
                          refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                          refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
                          refined_multipole.p_nonrot = &refined_multipole.p_rot;
                          refined_multipole.p_rot.axt = mpp.axt;
                          found_par_b = found_par_c = true;
                          refined_multipole.npole++;
                          break;
                        }
                        else
                        {
                          for (auto const d : coords.atoms(b).bonds())
                          {
                            if (d != c && d != a)
                            {
                              std::size_t const type_of_d(params.type(coords.atoms(d).energy_type(), tinker::MULTIPOLE));
                              if (type_of_d == mpp.index[3])
                              {
                                refined_multipole.axes.y() = d;
                                //refined_multipole.p_nonrot = &mpp; 
                                refined_multipole.p_rot.charge = mpp.charge;
                                refined_multipole.p_rot.dipole = mpp.dipole;
                                refined_multipole.p_rot.quadrupole[0U] = mpp.quadrupole[0U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[1U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[2U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[3U] = mpp.quadrupole[1U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[4U] = mpp.quadrupole[2U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[5U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[6U] = mpp.quadrupole[3U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[7U] = mpp.quadrupole[4U] * (B*B) / 3.0;
                                refined_multipole.p_rot.quadrupole[8U] = mpp.quadrupole[5U] * (B*B) / 3.0;
                                refined_multipole.p_nonrot = &refined_multipole.p_rot;
                                refined_multipole.p_rot.axt = mpp.axt;
                                found_par_b = found_par_c = true;
                                refined_multipole.npole++;
                                break;
                              }
                            }
                          }
                          if (found_par_c) break;
                        }
                      }
                    }
                  }
                }
              }
              if (found_par_b) break;
            }
          }
        }
        refined_multipole.p_rot.dipole *= B;

        if (refined_multipole.p_nonrot) 
          m_multipole.push_back(refined_multipole);
      }
    }


  }
  return m_multipole;
  //m_multipole_vec.push_back(m_multipole);
}

void tinker::refine::purge_nb_at_same_molecule(coords::Coordinates const & coords, tinker::refine::refined & params)
{
  
  auto& matrices = params.set_pair_matrices();
  for (std::size_t i = 0u; i < matrices.size(); i++)
  {
    for (std::size_t j = 0u; j < matrices.at(i).pair_matrix.rows(); j++)
    {
      for (std::size_t k = 0u; k < matrices.at(i).pair_matrix.cols(); k++)
      {
        auto& list = matrices.at(i).pair_matrix(j, k);
        std::vector<std::size_t> to_be_removed;
        for (std::size_t l = 0u; l < list.size(); l++)
        {
          auto& element = matrices.at(i).pair_matrix(j, k).at(l);
          bool found1 = false, found2 = false;
          std::size_t mol1(0), mol2(0);
          for (std::size_t m = 0u; m < coords.molecules().size(); m++)
          {
            if (found1 && found2)
              break;
            else
            {
              auto const& vec = coords.molecule(i);
              if (std::find(vec.begin(), vec.end(), element.a) != vec.end())
              {
                mol1 = m;
                found1 = true;
              }
              if (std::find(vec.begin(), vec.end(), element.b) != vec.end())
              {
                mol2 = m;
                found2 = true;
              }
              
            }
          }
          
          if (found1 && found2)
          {
            if (mol1 == mol2)
            {
              to_be_removed.push_back(l);
            }
          }
        }
        sort(to_be_removed.rbegin(), to_be_removed.rend());
        for (std::size_t l = 0u; l < to_be_removed.size(); l++)
        {
          list.erase(list.begin() + to_be_removed.at(l));
        }

      }
    }
  }


}

tinker::refine::vector_polarize tinker::refine::refine_pol(coords::Coordinates const & coords, tinker::parameter::parameters const & params)
{
  std::size_t const N(coords.atoms().size());
  refine::types::polarize refined_polarize;
  tinker::refine::vector_polarize m_polarize;
  for (std::size_t i(0U); i < N; i++)
  {
    for (auto const & pp : params.polarizes())
    {
      std::size_t const type_of_a(params.type(coords.atoms(i).energy_type(), tinker::POLARIZE));
      refined_polarize.atoms = pp.not_contracted_bond;

      if (type_of_a == pp.index)
      {
        refined_polarize.force = pp.p;
        refined_polarize.ff = pp.pp;

        refined_polarize.pdamp = pow(pp.p, (1.0 / 6.0));
        m_polarize.push_back(refined_polarize);
      }
    }
  }
  return m_polarize;
  //m_polarize_vec.push_back(m_polarize);
}

void tinker::refine::refined::clear(void)
{
  m_angles.clear();
  m_bonds.clear();
  m_impropers.clear();
  m_imptors.clear();
  m_torsions.clear();
  m_ureys.clear();
  m_multipole_vec.clear();
  m_polarize_vec.clear();
  m_pair_matrices.clear();
  for (auto & relation : m_relations) relation.clear();
  for (auto & vdwc_matrix : m_vdwc_matrices) vdwc_matrix.clear();
}

void tinker::refine::refined::swap_data(refined &rhs)
{
  m_angles.swap(rhs.m_angles);
  m_bonds.swap(rhs.m_bonds);
  m_impropers.swap(rhs.m_impropers);
  m_imptors.swap(rhs.m_imptors);
  m_torsions.swap(rhs.m_torsions);
  m_ureys.swap(rhs.m_ureys);
  m_pair_matrices.swap(rhs.m_pair_matrices);
  for (std::size_t i(0u); i < 5u; ++i) m_relations[i].swap(rhs.m_relations[i]);
  m_removes.swap(rhs.m_removes);
  m_red_types.swap(rhs.m_red_types);
  //m_multipole.swap(rhs.m_multipole);
  m_multipole_vec.swap(rhs.m_multipole_vec);
  m_polarize_vec.swap(rhs.m_polarize_vec);

  for (std::size_t i(0u); i < 6u; ++i) m_vdwc_matrices[i].swap(rhs.m_vdwc_matrices[i]);
}



std::ostream& tinker::refine::types::operator<< (std::ostream &stream, binary_quadratic const &bq)
{
  stream << "[" << bq.atoms[0] << "<->" << bq.atoms[1] << "]:(F: " << bq.force << ", L: " << bq.ideal << ")";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, ternary_quadratic const &bq)
{
  stream << "[" << bq.atoms[0] << "<->" << bq.atoms[1] << "<->" << bq.atoms[2] << "]:(F: " << bq.force << ", L: " << bq.ideal << ")";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, torsion const &bq)
{
  stream << "[" << bq.atoms[0] << "<->" << bq.atoms[1] << "<->" << bq.atoms[2] << "<->" << bq.atoms[3] << "]";
  for (std::size_t i(0u); i < bq.p.number; ++i)
  {
    stream << ":(F: " << bq.p.force[i] << ", A: " << bq.p.ideal[i] << ", O: " << bq.p.order[i] << ")";
  }
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, improper const &bq)
{
  stream << "[" << bq.center << "<.>" << bq.ligand[0] << "<.>" << bq.ligand[1] << "<.>" << bq.twist << "]";
  for (std::size_t i(0u); i < bq.p.number; ++i)
  {
    stream << ":(F: " << bq.p.force[i] << ", A: " << bq.p.ideal[i] << ", O: " << bq.p.order[i] << ")";
  }
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, imptor const &bq)
{
  stream << "[" << bq.center << "<.>" << bq.ligand[0] << "<.>" << bq.ligand[1] << "<.>" << bq.twist << "]";
  for (std::size_t i(0u); i < bq.p.number; ++i)
  {
    stream << ":(F: " << bq.p.force[i] << ", A: " << bq.p.ideal[i] << ", O: " << bq.p.order[i] << ")";
  }
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, multipole const &bq)
{
  stream << "MP:[" << bq.center << "<.>" << bq.axes.z() << "<.>" << bq.axes.x() << "<.>" << bq.axes.y() << "]\n";
  stream << ":(C: " << (*bq.p_nonrot).charge << ", DP: " << (*bq.p_nonrot).dipole << ")";
  stream << ":(QP: " << (*bq.p_nonrot).quadrupole(0) << " // " << (*bq.p_nonrot).quadrupole(1) << ", " << (*bq.p_nonrot).quadrupole(2) << " //";
  stream << "" << bq.p_rot.quadrupole(3) << ", " << bq.p_rot.quadrupole(4) << ", " << bq.p_rot.quadrupole(5) << ")\n";
  stream << ":(C: " << bq.p_rot.charge << ", DP: " << bq.p_rot.dipole << ")";
  stream << ":(QP: " << bq.p_rot.quadrupole(0) << " // " << bq.p_rot.quadrupole(1) << ", " << bq.p_rot.quadrupole(2) << " //";
  stream << "" << bq.p_rot.quadrupole(3) << ", " << bq.p_rot.quadrupole(4) << ", " << bq.p_rot.quadrupole(5) << ")";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, opbend const &bq)
{
  stream << "[" << bq.atoms[0] << "<->" << bq.atoms[1] << "<->" << bq.atoms[2] << "<->" << bq.atoms[3] << "]";
  stream << ":(F: " << bq.force << ")";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, polarize const &bq)
{
  stream << "[" << bq.center << "]";
  stream << ":(F: " << bq.force << ", F2: " << bq.ff << ")";
  stream << ":[" << bq.atoms[0] << "<->" << bq.atoms[1] << "<->" << bq.atoms[2] << "<->" << bq.atoms[3] << "]";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, strbend const &bq)
{
  stream << "[" << bq.atoms[0] << "<->" << bq.atoms[1] << "<->" << bq.atoms[2] << "<->" << bq.atoms[3] << "]";
  stream << ":(F: " << bq.force << ", F2: " << bq.ff << ")";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, nbpair const &bq)
{
  stream << "{" << std::setw(5) << bq.a << "," << std::setw(5) << bq.b << "}";
  return stream;
}

std::ostream& tinker::refine::types::operator<< (std::ostream &stream, nbpm const &bq)
{
  stream << "Pair Matrix with ID: " << bq.param_matrix_id << std::endl;
  for (std::size_t i(0u); i < bq.pair_matrix.size(); ++i)
  {
    for (auto pair : bq.pair_matrix(i))
    {
      stream << ", " << pair;
    }
    stream << std::endl;
  }
  return stream;
}

std::ostream& tinker::refine::operator<< (std::ostream &stream, refined const &ref)
{
  stream << "Bonds:" << std::endl;
  for (auto & a : ref.bonds()) stream << a << std::endl;
  stream << "Angles:" << std::endl;
  for (auto & a : ref.angles()) stream << a << std::endl;
  stream << "Impropers:" << std::endl;
  for (auto & a : ref.impropers()) stream << a << std::endl;
  stream << "Imptors:" << std::endl;
  for (auto & a : ref.imptors()) stream << a << std::endl;
  stream << "Torsions:" << std::endl;
  for (auto & a : ref.torsions()) stream << a << std::endl;
  stream << "Ureys:" << std::endl;
  for (auto & a : ref.ureys()) stream << a << std::endl;
  stream << "Pairs:" << std::endl;
  for (auto & a : ref.pair_matrices()) stream << a << std::endl;
  stream << "Multipoles:" << std::endl;
  for (auto mpv : ref.multipole_vecs())
  {
    for (auto & a : mpv) stream << a << std::endl;
  }
  stream << "Polarize:" << std::endl;
  for (auto pv : ref.polarize_vecs())
  {
    for (auto & a : pv) stream << a << std::endl;
  }
  stream << "VDWCM:" << std::endl;
  for (auto vdwcm : ref.vdwcm())
  {
    stream << "Matrix:" << std::endl;
    stream << vdwcm;
  }
  return stream;
}
