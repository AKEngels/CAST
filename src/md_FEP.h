#pragma once
#include "md.h"

namespace md
{
  /** collection of variables for FEP calculation
*/
  struct fepvar
  {
    /**lambda_el of former window for appearing atoms*/
    double mein;
    /**lambda_el of former window for disappearing atoms*/
    double meout;
    /**lambda_vdw of former window for appearing atoms*/
    double mvin;
    /**lambda_vdw of former window for disappearing atoms*/
    double mvout;
    /**lambda_el for appearing atoms*/
    double ein;
    /**lambda_el for disappearing atoms*/
    double eout;
    /**lambda_vdw for appearing atoms*/
    double vin;
    /**lambda_vdw for disappearing atoms*/
    double vout;
    /**lambda_el of next window for appearing atoms*/
    double dein;
    /**lambda_el of next window for disappearing atoms*/
    double deout;
    /**lambda_vdw of next window for appearing atoms*/
    double dvin;
    /**lambda_vdw of next window for disappearing atoms*/
    double dvout;
  };




}