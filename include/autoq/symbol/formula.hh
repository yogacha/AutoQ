#ifndef _AUTOQ_FORMULA_HH_
#define _AUTOQ_FORMULA_HH_

#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>
#include "autoq/complex/complex.hh"
#include <boost/functional/hash.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <utility>

struct VectorHash
{
    std::size_t operator()(const std::vector<int> &v) const
    {
        std::size_t hash = v.size();
        for (int i : v)
        {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

namespace AUTOQ
{
    namespace Symbol
    {
        struct Formula;
        using coefficient_map = std::unordered_map<std::vector<int>, Complex::Complex, VectorHash>;
    }
}

// Macro MODULAR_ARITHMETIC
#define MODULAR(x, m) ((m + (x % m)) % m)

// Formula symbol
struct AUTOQ::Symbol::Formula
{
private:
    bool internal = false;
    // bool is_leaf_v = true;  // redundant?
public:
    int period;            // omega^period = -1
    std::vector<int> vars; // x_1, x_2, ..., x_n
    std::unordered_map<std::vector<int>, Complex::Complex, VectorHash> formula;

    Formula(int p, const std::vector<int> &v, const std::unordered_map<std::vector<int>, Complex::Complex, VectorHash> &f) : period(p), vars(v), formula(f) {}
    Formula(int p, const std::vector<int> &v) : period(p), vars(v), formula({}) {}
    Formula() : period(1), vars({}), formula({}) {}
    Formula(Complex::Complex c) : period(1), vars({}), formula({{{}, c}}) {}

    static Formula Angle(const boost::rational<boost::multiprecision::cpp_int> &theta) { return Formula(Complex::Complex::Angle(theta)); }
    static Formula One() { return Formula(Complex::Complex::One()); }
    static Formula Zero() { return Formula(); }
    static Formula Rand() { return Formula(Complex::Complex::Rand()); }
    static Formula sqrt2() { return Formula(Complex::Complex::sqrt2()); }

    Formula to_period(int p) const
    {
        if (p % period != 0)
        {
            THROW_AUTOQ_ERROR("The period must be a multiple of the current period.");
        }
        Formula res(p, vars);
        int factor = p / period;
        std::vector<int> new_key(vars.size(), 0);

        for (auto &pair : formula)
        {
            for (size_t i = 0; i < vars.size(); i++)
            {
                new_key[i] = pair.first[i] * factor;
            }
            res.formula[new_key] = pair.second;
        }
        return res;
    }

    Formula to_variables(const std::vector<int> &v) const
    {
        Formula res(period, v);
        std::vector<int> index_map;
        std::vector<int> new_key(v.size(), 0);
        for (size_t i = 0; i < vars.size(); i++)
        {
            if (std::find(v.begin(), v.end(), vars[i]) != v.end())
            {
                index_map.push_back(std::distance(v.begin(), std::find(v.begin(), v.end(), vars[i])));
            }
            else
            {
                THROW_AUTOQ_ERROR("The new variables must be a subset of the current variables.");
            }
        }

        for (auto &[key, val] : formula)
        {
            for (size_t i = 0; i < key.size(); i++)
            {
                new_key[index_map[i]] = key[i];
            }
            res.formula[new_key] = val;
        }
        return res;
    }

    // unify two formulas, so that they have the same period and variables
    static std::pair<Formula, Formula> unify(const Formula &fa, const Formula &fb)
    {
        if (fa.vars != fb.vars)
        {
            std::vector<int> union_vars, va = fa.vars, vb = fb.vars;
            std::sort(va.begin(), va.end());
            std::sort(vb.begin(), vb.end());
            std::set_union(va.begin(), va.end(), vb.begin(), vb.end(), std::back_inserter(union_vars));
            return unify(fa.to_variables(union_vars), fb.to_variables(union_vars));
        }
        if (fa.period > fb.period)
        {
            return {fa, fb.to_period(fa.period)};
        }
        else if (fa.period < fb.period)
        {
            return {fa.to_period(fb.period), fb};
        }
        return {fa, fb};
    }

    Formula &fraction_simplification()
    {
        for (auto &[key, val] : formula)
        {
            formula[key] = val.fraction_simplification();
        }
        return *this;
    }
    Formula &key_simplification()
    {
        // erase zero coefficients and take modulo 2*period for the keys
        int mod = 2 * period;
        std::vector<int> new_key(vars.size(), 0);
        AUTOQ::Complex::Complex val;
        for (auto it = formula.begin(); it != formula.end();)
        {
            if (it->second.isZero())
            {
                it = formula.erase(it);
            }
            else if (std::all_of(it->first.begin(), it->first.end(), [mod](int i)
                                 { return 0 <= i && i < mod; }))
            {
                it++;
            }
            else
            {
                val = it->second;
                for (size_t i = 0; i < it->first.size(); i++)
                {
                    new_key[i] = MODULAR(it->first[i], mod);
                }
                it = formula.erase(it);
                val = val + formula[new_key];
                if (val.isZero())
                {
                    formula.erase(new_key);
                }
                else
                {
                    formula[new_key] = val;
                }
            }
            // TODO: do we need val.frac_simplification()?
        }
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const Formula &f)
    {
        os << "ω^" << 2 * f.period << " = 1; vars=";
        for (int i : f.vars)
        {
            os << "x" << i << ",";
        }
        os << std::endl;
        os << "formula= ";
        for (auto &[key, val] : f.formula)
        {
            os << val << " ω^";
            for (size_t i = 0; i < f.vars.size(); i++)
            {
                os << key[i] << "x" << f.vars[i] << "+";
            }
            os << "  ";
        }
        return os;
    }

    bool operator==(const Formula &o) const
    {
        if (period != o.period || vars != o.vars)
        {
            auto [fa, fb] = unify(*this, o);
            return fa == fb;
        }
        else
        {
            return Formula(*this).key_simplification().formula == Formula(o).key_simplification().formula;
        }
    }
    bool operator!=(const Formula &o) const { return !(*this == o); }
    // NOTE: simplicity of key is preserved under +, -, *, to_period, to_variables, negation, ite
    // so no need to call `key_simplification()` after these operations
    Formula operator+(const Formula &o) const
    {
        auto [fa, fb] = unify(*this, o);
        for (auto &[key, val] : fb.formula)
        {
            val = fa.formula[key] + val;
            if (val.isZero())
            {
                fa.formula.erase(key);
            }
            else
            {
                fa.formula[key] = val;
            }
        }
        return fa;
    }
    Formula operator-(const Formula &o) const
    {
        auto [fa, fb] = unify(*this, o);
        for (auto &[key, val] : fb.formula)
        {
            val = fa.formula[key] - val;
            if (val.isZero())
            {
                fa.formula.erase(key);
            }
            else
            {
                fa.formula[key] = val;
            }
        }
        return fa;
    }
    Formula operator*(const Formula &o) const
    {
        auto [fa, fb] = unify(*this, o);
        Formula res(fa.period, fa.vars);
        auto mod = 2 * fa.period;
        std::vector<int> new_key(fa.vars.size(), 0);
        for (auto &[key, val] : fa.formula)
        {
            for (auto &[key2, val2] : fb.formula)
            {
                for (size_t i = 0; i < key.size(); i++)
                {
                    new_key[i] = MODULAR(key[i] + key2[i], mod);
                }
                val = val * val2;
                if (val.isZero())
                {
                    res.formula.erase(new_key);
                }
                else
                {
                    res.formula[new_key] = val;
                }
            }
        }
        return res;
    }
    Formula negation(int x) const
    {
        auto it = std::find(vars.begin(), vars.end(), x);
        if (it == vars.end())
        {
            THROW_AUTOQ_ERROR("The variable x must be in the variables.");
        }
        auto n = std::distance(vars.begin(), it);
        int mod = 2 * period;
        std::vector<int> new_key;
        Formula res(period, vars);
        for (auto &[key, val] : formula)
        {
            new_key = key;
            new_key[n] = MODULAR(mod - key[n], mod);
            res.formula[new_key] = Complex::Complex(val).counterclockwise(
                boost::rational<boost::multiprecision::cpp_int>(key[n], mod));
        }
        return res;
    }
    static Formula ite(int x, const Formula &f1, const Formula &f2)
    {
        auto [fa, fb] = unify(f1, f2);
        if (std::find(fa.vars.begin(), fa.vars.end(), x) == fa.vars.end())
        {
            THROW_AUTOQ_ERROR("The variable x must be in the variables.");
        }
        return ( (fb + fa) + (fb - fa) * Formula(fa.period, {x}, {{{fa.period}, Complex::Complex(1)}})
            ).divide_by_the_square_root_of_two(2).fraction_simplification();
    }

    Formula &divide_by_the_square_root_of_two(int times = 1)
    {
        for (auto &[key, val] : formula)
        {
            formula[key] = val.divide_by_the_square_root_of_two(times);
        }
        return *this;
    }
    Formula &counterclockwise(const boost::rational<boost::multiprecision::cpp_int> &theta)
    {
        for (auto &[key, val] : formula)
        {
            formula[key] = val.counterclockwise(theta);
        }
        return *this;
    }
    Formula &clockwise(const boost::rational<boost::multiprecision::cpp_int> &theta)
    {
        for (auto &[key, val] : formula)
        {
            formula[key] = val.clockwise(theta);
        }
        return *this;
    }

    Formula evaluate(const std::vector<std::optional<bool>> &assign) const
    {
        if (assign.size() != vars.size())
        {
            THROW_AUTOQ_ERROR("The size of the assignment must be equal to the size of the variables.");
        }
        int k;
        std::vector<int> new_key;
        auto res = Formula(period, vars);
        for (auto &[key, val] : formula)
        {
            k = 0;
            new_key = key;
            for (size_t i = 0; i < assign.size(); i++)
            {
                if (assign[i].has_value())
                {
                    k += (assign[i].value() ? key[i] : 0);
                    new_key[i] = 0;
                }
            }
            res.formula[new_key] = res.formula[new_key] + val * Complex::Complex::Angle(
                boost::rational<boost::multiprecision::cpp_int>(k, 2 * period));
        }
        return res;
    }

    Formula qft(const std::vector<int> &qubits) const
    {
        std::vector<int> eval_var;
        for (size_t i = 0; i < qubits.size(); i++)
        {
            if (std::find(vars.begin(), vars.end(), qubits[i]) == vars.end())
            {
                THROW_AUTOQ_ERROR("The variables to be evaluated must be in the variables.");
            }
            else
            {
                eval_var.push_back(qubits[i]);
            }
        }

        Formula res(period, vars);
        return res;
    }

    bool isZero() const
    {
        if (formula.empty())
        {
            return true;
        }
        else
        {
            return Formula(*this).key_simplification().formula.empty();
        }
    }
};

#endif