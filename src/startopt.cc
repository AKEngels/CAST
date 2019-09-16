#include <memory>
#include "startopt.h"
#include "startopt_ringsearch.h"
#include "startopt_solvadd.h"

void startopt::apply(coords::Coordinates& c, coords::Ensemble_PES& e)
{
  std::unique_ptr<startopt::Preoptimizer> optimizer;
  std::size_t multi(Config::get().startopt.number_of_structures / e.size());
  std::cout << Config::get().startopt;
  std::cout << "-------------------------------------------------\n";
  if (Config::get().startopt.type == config::startopt::types::T::SOLVADD)
  {
    optimizer = std::unique_ptr<startopt::Preoptimizer>
    {
      new startopt::preoptimizers::Solvadd(c,
        Config::get().startopt.solvadd.maxDistance)
    };
  }
  else if (Config::get().startopt.type == config::startopt::types::T::RINGSEARCH ||
    Config::get().startopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
  {
    optimizer = std::unique_ptr<startopt::Preoptimizer>
    {
      new startopt::preoptimizers::R_evolution(c)
    };
    multi = Config::get().startopt.ringsearch.population;
  }
  if (optimizer != nullptr)
  {
    // Generate structures using initial preoptimizer (ringsearch or solvadd)
    //std::cout << "PreGenerate.\n";
    optimizer->generate(e, (multi > 0u ? multi : 1u));
    //std::cout << "PostGenerate.\n";
    if (Config::get().startopt.type == config::startopt::types::T::RINGSEARCH_SOLVADD)
    {
      std::size_t sa_multi(Config::get().startopt.number_of_structures / optimizer->PES().size());

      auto sap = std::unique_ptr<startopt::Preoptimizer>
      {
        new startopt::preoptimizers::Solvadd(optimizer->final_coords(),
          Config::get().startopt.solvadd.maxDistance)
      };
      sap->generate(optimizer->PES(), sa_multi);
      sap.swap(optimizer);
    }
    c.swap(optimizer->final_coords());
    e = optimizer->PES();
  }
}
