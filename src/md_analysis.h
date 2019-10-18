/**
CAST 3
md_analysis.h
Purpose: header for analysis of molecular dynamics simulation

@author Susanne Sauer
@version 1.0
*/

#pragma once 

#include <vector>
#include <string>

/**forward declaration of MD simulation class*/
namespace md {
  class simulation;
}

/**namespace for stuff that is used for analysis of MD simulation*/
namespace md_analysis
{
  /**struct that contains information about an atom pair that is to be analyzed*/
  struct ana_pair
  {
    /**index of first atom (starting with 0)*/
    int a;
    /**index of second atom (starting with 0)*/
    int b;
    /**element symbol of first atom*/
    std::string symbol_a;
    /**element symbol of second atom*/
    std::string symbol_b;
    /**name of first atom (element and tinker atom index)*/
    std::string name_a;
    /**name of second atom (element and tinker atom index)*/
    std::string name_b;
    /**legend of this atom pair in the graph*/
    std::string legend;
    /**distances for every MD frame*/
    std::vector<double> dists;

    /**constructor
    @param p1: tinker atom index of first atom (i.e. starting with 1)
    @param p2: tinker atom index of second atom (i.e. starting with 1)*/
    ana_pair(int p1, int p2) { a = p1 - 1; b = p2 - 1; }
  };

  /**information about a zone for which temperature is to by analyzed*/
  struct zone
  {
    /**legend for plotting*/
    std::string legend;
    /**atom indizes (starting with 0)*/
    std::vector<size_t> atoms;
    /**temperatures for every MD step*/
    std::vector<double> temperatures;
  };

  /**create the atom pairs to be analyzed*/
  void create_ana_pairs(md::simulation* md_obj);
  /**function that fills zones with atoms*/
  std::vector<zone> find_zones(md::simulation* md_obj);

  /**function to write distances into a file "distances.csv"
  @param pairs: atom pairs between which the distance should be calculated*/
  void write_dists_into_file(md::simulation* md_obj);
  /**function to write the temperatures for all zones into a file "zones.csv"*/
  void write_zones_into_file(md::simulation* md_obj);

  /**function to plot distances for atom pairs
    @param pairs: atom pairs to be plotted*/
  void plot_distances(md::simulation* md_obj);
  /**function to plot temperatures for all zones*/
  void plot_zones(md::simulation* md_obj);

  /**adding simulation information (run after every MD step)*/
  void add_analysis_info(md::simulation* md_obj);
  /**write information to file and plot it if possible (at the end of simulation)*/
  void write_and_plot_analysis_info(md::simulation* md_obj);
}