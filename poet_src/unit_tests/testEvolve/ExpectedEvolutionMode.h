/**\file
 *
 * \brief Declare & define a class tracking flags indicating the expected state
 * of the evolution (e.g. is wind saturated, current evolution mode).
 *
 * \ingroup UnitTests_group
 */

#ifndef __EXPECTED_EVOLUTION_MODE_H
#define __EXPECTED_EVOLUTION_MODE_H

#include "../shared/Common.h"
#include <list>

namespace Evolve {

    /**\brief Some evolution mode that changes at specified ages.
     *
     * \ingroup UnitTests_group
     */
    template<typename MODE_TYPE>
        class ExpectedEvolutionMode {
        private:
            ///The ages at which evolution mode changes occur.
            std::list<double> __age_breaks;

            ///The evolution modes that apply between consecutive __age_breaks.
            std::list<MODE_TYPE> __expected_mode;

            ///The precision with which breaks should be detected
            double __break_precision;
        public:
            ///Create.
            ExpectedEvolutionMode(double break_precision = 1e-5) :
                __break_precision(break_precision)
            {}

            ///Add an evolution mode that applies up to the given age.
            void add_break(double age, MODE_TYPE mode)
            {__age_breaks.push_back(age); __expected_mode.push_back(mode);}

            ///Is the given age close to a break (hence ambigous mode).
            bool near_break(double age) const
            {
                for(
                    std::list<double>::const_iterator 
                        break_i = __age_breaks.begin();
                    break_i != __age_breaks.end();
                    ++break_i
                )
                    if(check_diff(age, *break_i, 0.0, __break_precision))
                        return true;
                return false;
            }

            ///The evolution mode that corresponds to the given age.
            MODE_TYPE operator()(double age) const
            {
                typename std::list<MODE_TYPE>::const_iterator
                    mode_i = __expected_mode.begin();
                for(
                    std::list<double>::const_iterator age_i = __age_breaks.begin();
                        age_i != __age_breaks.end();
                        age_i++
                ) {
                    std::list<double>::const_iterator next_age_i=age_i;
                    ++next_age_i;
                    if(
                        *age_i <= age
                        &&
                        (next_age_i == __age_breaks.end() || *next_age_i > age)
                    )
                        return *mode_i;
                    ++mode_i;
                }
                std::ostringstream msg;
                msg << "Age "
                    << age
                    << " outside the range for which expected evolution modes "
                    << "are defined.";
                throw Core::Error::BadFunctionArguments(msg.str());
            }
        };



}//End Evolve namespace.

#endif
