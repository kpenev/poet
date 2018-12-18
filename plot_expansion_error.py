#!/usr/bin/env python3

"""Create plots showing the error due to truncating the ecc. expansion."""

from matplotlib import pyplot
import numpy

def plot_actual_error(data):
    """Create plots of the actual (fine-rough) error in rates."""

    for e_order in numpy.unique(data['e_order']):
        this_e_order = (data['e_order'] == e_order)
        pyplot.subplot(221)
        pyplot.semilogy(
            data['e'][this_e_order],
            numpy.abs(
                (
                    data['fine_a_rate'][this_e_order]
                    -
                    data['rough_a_rate'][this_e_order]
                )
                /
                data['fine_a_rate'][this_e_order]
            ),
            label='$\\mathcal{O}(e^{%d})$' % e_order
        )
        pyplot.title('$\\left|\\delta \\dot{a} / a\\right|$')
        pyplot.ylim(1e-7, 10)

        pyplot.subplot(222)
        pyplot.semilogy(
            data['e'][this_e_order],
            numpy.abs(
                (
                    data['fine_e_rate'][this_e_order]
                    -
                    data['rough_e_rate'][this_e_order]
                )
                /
                data['fine_e_rate'][this_e_order]
            )

        )
        pyplot.title('$\\left|\\delta \\dot{e} / e\\right|$')
        pyplot.ylim(1e-7, 10)

        pyplot.subplot(223)
        pyplot.semilogy(
            data['e'][this_e_order],
            numpy.abs(
                (
                    data['fine_spin_rate'][this_e_order]
                    -
                    data['rough_spin_rate'][this_e_order]
                )
                /
                data['fine_spin_rate'][this_e_order]
            )
        )
        pyplot.title(
            '$\\left|\\delta \\dot{\\Omega_\\star} / \\Omega_\\star\\right|$'
        )
        pyplot.ylim(1e-7, 10)

    pyplot.figlegend()

def plot_actual_vs_estimated_error(data, quantity, ratio=False):
    """Create plots comparing the actual and estimated errors."""

    e_orders = numpy.unique(data['e_order'])
    plot_columns = int(numpy.ceil(numpy.sqrt(e_orders.size)))
    plot_rows = int(numpy.ceil(e_orders.size / plot_columns))

    highest_e_order = (data['e_order'] == e_orders[-1])

    for plot_index, e_order in enumerate(e_orders):
        this_e_order = (data['e_order'] == e_order)
        pyplot.subplot(plot_rows, plot_columns, plot_index + 1)
        if ratio:
            pyplot.semilogy(data['e'][highest_e_order],
                            numpy.abs(
                                (
                                    data['fine_%s_rate' % quantity][highest_e_order]
                                    -
                                    data['rough_%s_rate' % quantity][this_e_order]
                                )
                                /
                                data['%s_rate_error' % quantity][this_e_order]
                            )
                            ,
                            label='$\\delta %s$: actual / estimated' % quantity)
            pyplot.ylim(0, 10)
        else:
            pyplot.semilogy(data['e'][this_e_order],
                            numpy.abs(
                                (
                                    data['fine_%s_rate' % quantity][this_e_order]
                                    -
                                    data['rough_%s_rate' % quantity][this_e_order]
                                )
                                /
                                data['fine_%s_rate' % quantity][this_e_order]
                            ),
                            label='fine - rough')

            pyplot.semilogy(data['e'][highest_e_order],
                            numpy.abs(
                                (
                                    data['fine_%s_rate' % quantity][highest_e_order]
                                    -
                                    data['rough_%s_rate' % quantity][this_e_order]
                                )
                                /
                                data['fine_%s_rate' % quantity][this_e_order]
                            ),
                            label='best - rough')

            pyplot.semilogy(data['e'][this_e_order],
                            (
                                data['%s_rate_error' % quantity][this_e_order]
                                /
                                numpy.abs(
                                    data['fine_%s_rate' % quantity][highest_e_order]
                                )
                            ),
                            label='estimated')
            pyplot.axhline(y=0.1)
            pyplot.ylim(1e-7, 10)
        pyplot.title('$\\mathcal{O}(e^{%d})$' % e_order)
    pyplot.figlegend(markerscale=5)

if __name__ == '__main__':
    plot_data = numpy.genfromtxt('error_test.txt', names=True, dtype=None)

    plot_actual_vs_estimated_error(plot_data, 'a', False)
    pyplot.show()

    plot_actual_vs_estimated_error(plot_data, 'e', False)
    pyplot.show()

    plot_actual_vs_estimated_error(plot_data, 'spin', False)
    pyplot.show()
