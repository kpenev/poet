/**\file
 *
 * \brief Declares a DissipatingBody with exponentially decaying differential
 * rotation coupling.
 *
 * \ingroup Star_group
 */


#ifndef __EXPONENTIAL_DECAY_DIFF_ROT_COUPLING_BODY_H
#define __EXPONENTIAL_DECAY_DIFF_ROT_COUPLING_BODY_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Evolve/DissipatingBody.h"
#include "../Evolve/ZoneOrientation.h"

namespace Star {

    ///\brief A body with differential rotation torque between two zones given
    ///by: \f$\dot{L}_1=\frac{I_1 L_2 - I_2 L_1}{\tau_c(I_1+I_2)}\f$.
    class LIB_PUBLIC ExponentialDecayDiffRotBody :
        virtual public Evolve::DissipatingBody {
    private:
        ///The age for which the last conigure() was called.
        double __current_age,
               
               ///The timescale for the differential rotation torque
               __timescale;

        ///\brief The coupling torque and its nonzero entries.
        ///
        ///The outer index is the index of the outer zone being torqued.
        ///The meaning of the inner index is defined by torque_entry.
        mutable std::valarray< std::valarray<Eigen::Vector3d> > __torque;

        ///Resets __torque to all undefined (correcting the size if necessary).
        void reset_torque();

        ///\brief Returns the entry in __torque that corresponds to the given
        ///quantity entry.
        ///
        ///The quantity to differentiate against must have an entry in __torque
        ///(i.e. the entry must not be zero in general).
        ///
        ///See DissipatingBody::angular_momentum_coupling() for a description of
        ///the arguments.
        Eigen::Vector3d &torque_entry(
            unsigned top_zone_index, 
            Evolve::Dissipation::QuantityEntry entry,
            bool wih_respect_to_top
        ) const;
    public:
        ///Construct a body with default differential rotatino torque.
        ExponentialDecayDiffRotBody(
                ///The timescale for differential rotation coupling.
                double coupling_timescale) :
            __current_age(Core::NaN),
            __timescale(coupling_timescale),
            __torque(0)
        {}

        ///\brief See DissipatingBody::configure().
        ///
        ///Invokes DissipatingBody::configure().
        virtual void configure(bool initialize,
                               double age,
                               double companion_mass,
                               double semimajor,
                               double eccentricity,
                               const double *spin_angmom,
                               const double *inclination = NULL,
                               const double *periapsis = NULL,
                               bool locked_surface = false,
                               bool zero_outer_inclination = false,
                               bool zero_outer_periapsis = false);

        ///See DissipatingBody::angular_momentum_coupling().
        Eigen::Vector3d angular_momentum_coupling(
            unsigned top_zone_index,
            Evolve::Dissipation::QuantityEntry
                entry=Evolve::Dissipation::NO_DERIV,
            bool with_respect_to_top=false
        ) const;
    };//End ExponentialDecayDiffRotBody class.

}//End Star namespace.

#endif
