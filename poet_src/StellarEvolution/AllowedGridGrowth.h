/**\file
 *
 * \brief Declares & defines a class describing how the mass-metallicity
 * interpolation grid is allowed to grow.
 */

#ifndef __ALLOWED_GRID_GROWTH_H
#define __ALLOWED_GRID_GROWTH_H

namespace StellarEvolution {

    ///The four directions in which an interpolation grid can grow.
    class AllowedGridGrowth {
    private:
        bool 
            ///Can grow toward lower masses?
            __lighter,

            ///Can grow toward higher masses?
            __heavier,

            ///Can grow toward lower metallicities?
            __poorer,

            ///Can grow toward higher metallicities?
            __richer;

    public:
        AllowedGridGrowth() :
            __lighter(true), __heavier(true), __poorer(true), __richer(true)
        {}

        ///Can grow toward lower masses?
        bool lighter() const {return __lighter;}

        ///Can grow toward higher masses?
        bool heavier() const {return __heavier;}

        ///Can grow toward lower metallicities?
        bool poorer() const {return __poorer;}

        ///Can grow toward higher metallicities?
        bool richer() const {return __richer;}

        ///Disable growth to lower masses.
        AllowedGridGrowth &block_lighter() {__lighter = false; return *this;}

        ///Disable growth to higher masses.
        AllowedGridGrowth &block_heavier() {__heavier = false; return *this;}

        ///Disable growth to lower metallicities.
        AllowedGridGrowth &block_poorer() {__poorer = false; return *this;}

        ///Disable growth to higher metallicities.
        AllowedGridGrowth &block_richer() {__richer = false; return *this;}

        ///Is growth allowed in at least one direction?
        operator bool() const
        {return __lighter || __heavier || __poorer || __richer;}
    }; //End AllowedGridGrowth class declaration.

} //End StellarEvolution namespace.

#endif
