// include/BoostIncludes.h
#ifndef BOOST_INCLUDES_H
#define BOOST_INCLUDES_H

// Boost math libraries
#include <boost/math/special_functions/bessel.hpp>                 // BOOST LIBRARIES:  1. BesselK in external fields
#include <boost/math/special_functions/beta.hpp>                   // Beta function for recursive relation in scatterred fields
#include <boost/math/special_functions/legendre.hpp>               // Lengendre Plm
#include <boost/math/quadrature/gauss_kronrod.hpp> 
#include <boost/math/quadrature/exp_sinh.hpp>

double Pi = boost::math::constants::pi<double>();

#endif