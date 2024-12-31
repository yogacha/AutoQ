/*****************************************************************************
 *  AUTOQ Tree Automata Library
 *
 *  Copyright (c) 2011  Ondra Lengal <ilengal@fit.vutbr.cz>
 *
 *  Description:
 *    Test suite for explicit tree automaton
 *
 *****************************************************************************/

// AUTOQ headers
#include <cmath>
#include <fstream>
#include <filesystem>

#include "autoq/error.hh"
#include "autoq/inclusion.hh"
#include "autoq/util/util.hh"
#include "autoq/complex/complex.hh"
#include "autoq/symbol/concrete.hh"
#include "autoq/aut_description.hh"
#include "autoq/parsing/timbuk_parser.hh"
#include "autoq/serialization/timbuk_serializer.hh"
#include "autoq/symbol/formula.hh"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AutType
#include <boost/test/unit_test.hpp>

using AUTOQ::Complex::Complex;
using AUTOQ::Symbol::coefficient_map;
using AUTOQ::Symbol::Concrete;
using AUTOQ::Symbol::Formula;

BOOST_AUTO_TEST_CASE(formula_variable_eq)
{
    Formula f1(4, {1, 2}, {{{1, 2}, AUTOQ::Complex::Complex::One()}});
    Formula f2(4, {2, 1, 3}, {{{2, 1, 0}, AUTOQ::Complex::Complex::One()}});
    BOOST_REQUIRE(f1.to_variables(f2.vars) == f2);
    BOOST_REQUIRE(f1 == f2);
}

BOOST_AUTO_TEST_CASE(formula_period_eq)
{
    Formula f1(4, {1, 2}, {{{1, 3}, AUTOQ::Complex::Complex::One()}});
    Formula f2(8, {1, 2}, {{{2, 6}, AUTOQ::Complex::Complex::One()}});
    BOOST_REQUIRE(f1.to_period(f2.period) == f2);
    BOOST_REQUIRE(f1 == f2);
}

BOOST_AUTO_TEST_CASE(formula_add)
{
    Formula f1(2, {1, 2}, {{{1, 2}, AUTOQ::Complex::Complex::One()}});
    Formula f2(4, {2, 1}, {{{4, 2}, AUTOQ::Complex::Complex::One()}});
    Formula f3(4, {1, 2}, {{{2, 4}, AUTOQ::Complex::Complex::One() + AUTOQ::Complex::Complex::One()}});
    BOOST_REQUIRE(f1 == f2);
    BOOST_REQUIRE(f1 + f2 == f3);
    auto f4 = f1 - f2;
}

BOOST_AUTO_TEST_CASE(formula_simplification)
{
    Formula f1(4, {1}, {
                           {{-5}, AUTOQ::Complex::Complex::One()}, // 3
                           {{21}, AUTOQ::Complex::Complex::One()}, // 5
                           {{11}, AUTOQ::Complex::Complex::One()}, // 3
                           {{4}, AUTOQ::Complex::Complex::Zero()}  // 4, drop
                       });
    BOOST_REQUIRE(f1.key_simplification() ==
                  Formula(4, {1}, {{{3}, AUTOQ::Complex::Complex(2)}, {{5}, AUTOQ::Complex::Complex::One()}}));
}

BOOST_AUTO_TEST_CASE(formula_mul)
{
    Formula f1(4, {1, 2}, {{{7, 2}, AUTOQ::Complex::Complex(3)}});
    Formula f2(2, {2, 1}, {{{3, 2}, AUTOQ::Complex::Complex(7)}});
    Formula f3(4, {1, 2}, {{{3, 0}, AUTOQ::Complex::Complex(21)}, {{1, 4}, AUTOQ::Complex::Complex(0)}});
    BOOST_REQUIRE(f1 * f2 == f3);
}

BOOST_AUTO_TEST_CASE(formula_negation)
{
    auto theta = boost::rational<boost::multiprecision::cpp_int>(3, 4);
    Formula f1(2, {1, 2}, {{{3, 2}, AUTOQ::Complex::Complex(1)}});
    Formula f2(2, {1, 2}, {{{1, 2}, AUTOQ::Complex::Complex::Angle(theta)}});
    BOOST_REQUIRE(f1.negation(1) == f2);
}

BOOST_AUTO_TEST_CASE(formula_ite)
{
    Formula f1 = Formula::ite(1,
                              Formula(2, {1}, {{{1}, AUTOQ::Complex::Complex(-1)}}),
                              Formula(2, {1}, {{{2}, AUTOQ::Complex::Complex::One()}}));
    Formula f2 = Formula(2, {1}, {{{3}, AUTOQ::Complex::Complex::One()}});

    BOOST_REQUIRE(f1.evaluate({true}).fraction_simplification() == f2.evaluate({true}));
    BOOST_REQUIRE(f1.evaluate({false}).fraction_simplification() == f2.evaluate({false}));
}

BOOST_AUTO_TEST_CASE(formula_qft)
{
    Formula f1 = Formula::One().to_variables({1, 3});
    Formula f2 = f1.to_variables({1, 2, 3}).qft({3, 1}).to_period(8);
    Formula f3 = f1.to_period(8).qft({3, 1});
    std::cout << f2 << std::endl;
    std::cout << f3 << std::endl;
    BOOST_REQUIRE(f2 == f3);
}