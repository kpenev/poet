#!/usr/bin/env python3

"""Utilities for creating custom MIST tracks on the fly."""

from matplotlib import pyplot
from mesa_reader import MesaLogDir
from astropy import units

def process_profile(profile_data):
    """Return tachocline M and R, and moments of inertia for given profile."""

    conv_bottom_i = 0
    while profile_data.mixing_type[conv_bottom_i] > 0:
        conv_bottom_i += 1
    dinertia = (2.0 / 3.0
                *
                profile_data.rmid**2
                *
                (profile_data.dm * units.g).to_value(units.M_sun))
    return (
        profile_data.mass[conv_bottom_i],
        profile_data.radius[conv_bottom_i],
        dinertia[:conv_bottom_i].sum(),
        dinertia[conv_bottom_i:].sum(),
    )


if __name__ == '__main__':
    mesa_data = MesaLogDir(
        '/home/kpenev/projects/MESA/mesa-r7503/star/test_mist/LOGS'
    )
    history = mesa_data.history_data
    pyplot.semilogx(history.star_age, history.log_R)
    pyplot.title('$R_\star(t)$')
    pyplot.show()
    pyplot.clf()

    for model_number in mesa_data.model_numbers:
        profile_age = mesa_data.history_data.star_age[model_number - 1]
        profile = mesa_data.profile_data(model_number=model_number)
    #    print('Mcore=%.3f, Rcore=%.3f, Ienv=%.3e, Icore=%.3e'
    #          %
    #          process_profile(profile))
        pyplot.plot(
            profile.radius / 10.0**history.log_R[model_number - 1],
            profile.mixing_type,
            '-g'
        )
        pyplot.axvline(
            history.rcore[model_number - 1] / 10.0**history.log_R[model_number - 1],
            color='green'
        )
        pyplot.plot(profile.mass, profile.mixing_type, '-r')
        pyplot.axvline(history.mcore[model_number - 1], color='red')
        pyplot.title('Mixing($M$ and $R$) at $t=%.3e$ Gyr' % (profile_age/1e9))
        pyplot.show()
        pyplot.clf()
