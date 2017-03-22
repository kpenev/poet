#!/usr/bin/python3

import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot

import sys
sys.path.append('../PythonPackage')

from interpolator_manager_gui import InterpolatorManagerGUI
from stellar_evolution.library_interface import\
    library as stellar_evolution_library
from matplotlib.backends.backend_tkagg import\
    FigureCanvasTkAgg,\
    NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from glob import glob
import os.path
import re
import numpy
import astropy
from datetime import datetime, timedelta
import functools
from math import floor

if sys.version_info[0] < 3:
    import Tkinter as Tk
    from Tkinter import ttk
else:
    import tkinter as Tk
    from tkinter import ttk

serialized_interpolator_dir = '../stellar_evolution_interpolators'
mesa_track_dir = '../MESA_tracks'

class InterpolatedQuantitySum :

    def __init__(self, q1, q2) :
        self.q1 = q1
        self.q2 = q2
        self.min_age = max(q1.min_age, q2.min_age)
        self.max_age = min(q1.max_age, q2.max_age)

    def __call__(self, age) :
        return self.q1(age) + self.q2(age)

    def deriv(self, age) :
        return self.q1.deriv(age) + self.q2.deriv(age)

class InterpolationInteractive :

    def quit(self) :
        """Exit the application."""

        self.window.quit()     # stops mainloop
        self.window.destroy()  # this is necessary on Windows to prevent
                               # Fatal Python Error: PyEval_RestoreThread:
                               # NULL tstate

    def get_interpolated_quantity(self, star_mass, star_metallicity) :
        """Return a callable for plotting an interpolation quantity."""

        if (
                star_mass < self.track_mass[0]
                or
                star_mass > self.track_mass[-1]
                or
                star_metallicity < self.track_metallicity[0]
                or
                star_metallicity > self.track_metallicity[-1]
        ) :
            return None

        if self.plot_quantity == 'I' :
            return InterpolatedQuantitySum(
                self.interpolator_manager.current_interpolator()(
                    'ICONV',
                    star_mass,
                    star_metallicity
                ),
                self.interpolator_manager.current_interpolator()(
                    'IRAD',
                    star_mass,
                    star_metallicity
                )
            )
        else :
            if self.plot_quantity == 'R' : quantity_id = 'radius'
            else : quantity_id = self.plot_quantity
            return self.interpolator_manager.current_interpolator()(
                quantity_id,
                star_mass,
                star_metallicity
            )

    def plot_interpolation(self,
                           star_mass,
                           star_metallicity,
                           plot_func,
                           deriv_order,
                           single_plot = False) :
        """Plot an interpolated stellar evolution track."""

        interpolated_quantity = self.get_interpolated_quantity(
            star_mass,
            star_metallicity
        )

        try : resolution =  int(self.resolution_text.get())
        except : resolution = 100
        interpolation_ages = numpy.exp(
            numpy.linspace(
                numpy.log(max(interpolated_quantity.min_age, 1e-5)),
                numpy.log(interpolated_quantity.max_age),
                resolution
            )[1:-1]
        )

        plot_x = age_transform(star_mass,
                               star_metallicity,
                               interpolation_ages)
        if deriv_order > 0 :
            derivatives = interpolated_quantity.deriv(interpolation_ages)
            plot_y = derivatives[deriv_order if single_plot else 0]
        else :
            plot_y = interpolated_quantity(interpolation_ages)

        plot_func(plot_x, plot_y, self.interp_plot_style.get())

        if single_plot : return

        if deriv_order > 0 :
            d1_y = numpy.copy(derivatives[1])
            if deriv_order > 1 :
                d2_y = numpy.copy(derivatives[2])

            if self.logy : 
                d1_y /= plot_y
                if deriv_order > 1 :
                    d2_y = (d2_y / plot_y - d1_y**2)

            if self.logx : 
                deriv_plot_funcname = 'semilogx'
                d1_y *= plot_x
                if deriv_order > 1 :
                    d2_y = d1_y + plot_x**2 * d2_y
            else :
                deriv_plot_funcname = 'plot'


            getattr(self.first_deriv_axes, deriv_plot_funcname)(
                plot_x,
                d1_y,
                self.interp_plot_style.get()
            )
            if deriv_order > 1 :
                getattr(self.second_deriv_axes, deriv_plot_funcname)(
                    plot_x,
                    d2_y,
                    self.interp_plot_style.get()
                )

    def separate_window(self, deriv_order) :
        """Spawn a new window plotting curves of given order derivative."""

        if self.logx and self.logy : plot = matplotlib.pyplot.loglog
        elif self.logx and not self.logy : plot = matplotlib.pyplot.semilogx
        elif not self.logx and self.logy : plot = matplotlib.pyplot.semilogy
        else : plot = matplotlib.pyplot.plot
        self.plot_interpolation(self.interp['mass'],
                                self.interp['metallicity'],
                                plot,
                                deriv_order,
                                True)
        matplotlib.pyplot.show()

    def display(self) :
        """(Re-)draw the plot as currently configured by the user."""

        if self.do_not_display : return

        exec('def age_transform(m, feh, t) : return '
             +
             self.age_transform_entry.get(),
             globals())
        main_x_lim = self.main_axes.get_xlim()
        main_y_lim = self.main_axes.get_ylim()

        self.main_axes.cla()
        if self.first_deriv_axes is not None :
            first_deriv_xlim = self.first_deriv_axes.get_xlim()
            first_deriv_ylim = self.first_deriv_axes.get_ylim()
            self.first_deriv_axes.cla()
        if self.second_deriv_axes is not None :
            second_deriv_xlim = self.second_deriv_axes.get_xlim()
            second_deriv_ylim = self.second_deriv_axes.get_ylim()
            self.second_deriv_axes.cla()

        if self.logx and self.logy : plot = self.main_axes.loglog
        elif self.logx and not self.logy : plot = self.main_axes.semilogx
        elif not self.logx and self.logy : plot = self.main_axes.semilogy
        else : plot = self.main_axes.plot

        if self.interpolate :
            self.plot_interpolation(self.interp['mass'],
                                    self.interp['metallicity'],
                                    plot,
                                    self.deriv)

        for track_mass_index, track_mass in enumerate(self.track_mass) :
            nearby_track_mass = (
                track_mass_index == self.track_below['mass']
                or
                track_mass_index == self.track_below['mass'] + 1
            )
            for track_metallicity_index, track_metallicity in enumerate(
                    self.track_metallicity
            ) :
                nearby_track_metallicity = (
                    track_metallicity_index == (
                        self.track_below['metallicity']
                    )
                    or
                    track_metallicity_index == (
                        self.track_below['metallicity'] + 1
                    )
                )
                if (
                        self.track_state[track_mass][track_metallicity].get()
                        or
                        (nearby_track_mass and nearby_track_metallicity)
                ) :
                    track_mass = self.track_mass[track_mass_index]
                    track_metallicity = self.track_metallicity[
                        track_metallicity_index
                    ]

                    if (
                            self.track_state[track_mass][track_metallicity]
                            and
                            self.interpolate
                    ) :
                        self.plot_interpolation(track_mass,
                                                track_metallicity,
                                                plot,
                                                self.deriv)

                    plot(age_transform(track_mass,
                                       track_metallicity,
                                       self.tracks[
                                           track_mass
                                       ][
                                           track_metallicity
                                       ][
                                           't'
                                       ]),
                         self.tracks[
                             track_mass
                         ][
                             track_metallicity
                         ][
                             self.plot_quantity
                         ],
                         self.track_plot_style.get())
        if not self.main_auto_axes :
            self.main_axes.set_xlim(main_x_lim)
            self.main_axes.set_ylim(main_y_lim)
        if self.first_deriv_axes is not None :
            if self.first_deriv_auto_axes :
                self.first_deriv_axes.set_xlim(self.main_axes.get_xlim())
            else :
                self.first_deriv_axes.set_xlim(first_deriv_xlim)
                self.first_deriv_axes.set_ylim(first_deriv_ylim)
            if self.second_deriv_axes is not None :
                if self.second_deriv_auto_axes :
                    self.second_deriv_axes.set_xlim(
                        self.main_axes.get_xlim()
                    )
                else :
                    self.second_deriv_axes.set_xlim(second_deriv_xlim)
                    self.second_deriv_axes.set_ylim(second_deriv_ylim)
        self.main_auto_axes = False
        self.first_deriv_auto_axes = False
        self.second_deriv_auto_axes = False
        self.main_canvas.show()

    def on_key_event(self, event) :
        """Handle standard matplotlib key presses."""

        key_press_handler(event, self.main_canvas, self.toolbar)

    def toggle_log_y(self) :
        """Switch between log and linear scale for the y axis and refresh."""

        self.logy = not self.logy
        self.logy_button.config(
            relief = Tk.SUNKEN if self.logy else Tk.RAISED
        )
        self.main_auto_axes = True
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.display()

    def toggle_log_x(self) :
        """Switch between log and linear scale for the x axis and refresh."""

        self.logx = not self.logx
        self.logx_button.config(
            relief = Tk.SUNKEN if self.logx else Tk.RAISED
        )
        self.main_auto_axes = True
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.display()

    def toggle_interpolation(self) :
        """Switch between displaying and not an interpolated track."""

        self.interpolate = not self.interpolate
        self.interpolate_button.config(
            relief = Tk.SUNKEN if self.interpolate else Tk.RAISED
        )
        self.display()

    def toggle_deriv(self, order) :
        """React appropriately toggling the derivative of the given order."""

        if self.deriv >= order : self.deriv = order - 1
        else : self.deriv = order

        if self.first_deriv_axes is not None :
            self.main_figure.delaxes(self.first_deriv_axes)
            self.first_deriv_axes = None
        if self.second_deriv_axes is not None :
            self.main_figure.delaxes(self.second_deriv_axes)
            self.second_deriv_axes = None

        if self.deriv == 0 :
            self.main_axes.set_position([0.05, 0.05, 0.9, 0.9])
        else :
            self.main_axes.set_position([0.05, 0.5, 0.9, 0.45])
            if self.deriv == 1 :
                self.first_deriv_axes = self.main_figure.add_axes(
                    [0.05, 0.05, 0.9, 0.45]
                )
            else :
                assert(self.deriv == 2)
                self.first_deriv_axes = self.main_figure.add_axes(
                    [0.05, 0.05, 0.45, 0.45]
                )
                self.second_deriv_axes = self.main_figure.add_axes(
                    [0.5, 0.05, 0.45, 0.45]
                )

        self.first_deriv_button.config(
            relief = Tk.RAISED if self.deriv < 1 else Tk.SUNKEN
        )
        self.second_deriv_button.config(
            relief = Tk.RAISED if self.deriv < 2 else Tk.SUNKEN
        )
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.display()

    def change_plot_quantity(self, plot_quantity = None) :
        """Modify y axis controls for the new quantity and refresh."""

        if plot_quantity is None :
            plot_quantity = self.selected_plot_quantity.get()
        self.plot_quantity = plot_quantity

        self.first_deriv_button.config(
            state = (Tk.DISABLED if self.max_deriv[plot_quantity] == 0
                     else Tk.NORMAL)
        )
        self.second_deriv_button.config(
            state = (Tk.DISABLED if self.max_deriv[plot_quantity] < 2
                     else Tk.NORMAL)
        )
        self.main_auto_axes = True
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.change_interpolated_quantity()
        self.display()

    def change_interpolated_quantity(self) :
        """Re-creates the interpolated quantity per the current settings."""


    def change_interp(self, quantity, input_new_value) :
        """
        Modify the stellar mass or [Fe/H] for the displayed interpolation.

        Args:
            - quantity: One of 'mass' or 'metallicity'.
            - input_new_value: The new value to set. Should be convertible
                               to float.

        Returns: 
            - None: if the quantity was not changed due to it still being
                    within the last change timeout.
            - True: if the quantity was actually changed.
        """

        if self.changing[quantity] : return

        new_value = float(input_new_value)
        self.interp[quantity] = new_value
        track_quantities = getattr(self, 'track_' + quantity)
        index_below = len(track_quantities) - 2
        while (index_below >= 0
               and
               track_quantities[index_below] > new_value
        ) : index_below -= 1
        self.track_below[quantity] = index_below
        self.track_below[quantity] = -2

        self.fine_interp_scale[quantity].config(
            from_ = (track_quantities[index_below] + 1e-5),
            to = (track_quantities[index_below + 1] - 1e-5)
        )

        self.coarse_interp_scale[quantity].set(new_value)
        self.fine_interp_scale[quantity].set(new_value)

        self.change_interpolated_quantity()
        self.changing[quantity] = True

        if self.display_job : self.window.after_cancel(self.display_job)
        self.display_job = self.window.after(10, self.display)

        self.changing[quantity] = False
        return True

    def toggle_all_tracks(self) :
        """Enable displaying all tracks."""

        selected = (self.all_tracks_button.cget('relief') == Tk.RAISED)
        self.all_tracks_button.config(
            relief = (Tk.SUNKEN if selected else Tk.RAISED),
            text = ('None' if selected else 'All')
        )

        for button in self.track_mass_button.values() :
            button.config(relief = (Tk.SUNKEN if selected else Tk.RAISED))
        for button in self.track_metallicity_button.values() :
            button.config(relief = (Tk.SUNKEN if selected else Tk.RAISED))

        for mass_track_state in self.track_state.values() :
            for track_state in mass_track_state.values() :
                track_state.set(selected)

        self.display()

    def toggle_track_mass(self, mass) :
        """Select/deselect displaying all tracks with the given mass."""

        selected = (self.track_mass_button[mass].cget('relief') == Tk.RAISED)
        self.track_mass_button[mass].config(
            relief = (Tk.SUNKEN if selected else Tk.RAISED)
        )
        for track_state in self.track_state[mass].values() :
            track_state.set(selected)

        self.display()

    def toggle_track_metallicity(self, metallicity) :
        """Select/deselect displaying all tracks with the given [Fe/H]."""

        selected = (
            self.track_metallicity_button[metallicity].cget('relief')
            == 
            Tk.RAISED
        )
        self.track_metallicity_button[metallicity].config(
            relief = (Tk.SUNKEN if selected else Tk.RAISED)
        )
        for mass_track_state in self.track_state.values() :
            if metallicity in mass_track_state :
                mass_track_state[metallicity].set(selected)

        self.display()

    def auto_axes(self) :
        """Re-plot letting matplotlib determine the axes limits."""

        self.main_auto_axes = True
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.display()

    def __init__(self, window, tracks_dir) :
        """Setup user controls and display frame."""

        def create_main_axes() :
            """Create a figure and add an axes to it for drawing."""

            self.main_figure = Figure(figsize=(5, 4), dpi=100)
            self.main_axes = self.main_figure.add_axes(
                (0.05, 0.05, 0.9, 0.9)
            )
            self.first_deriv_axes = None
            self.second_deriv_axes = None

        def create_axes_controls() :
            """Create controls for plot quantity and log axes."""

            def create_x_controls() :
                """Create the controls for the x axis."""

                x_controls_frame = Tk.Frame(window)
                x_controls_frame.grid(row = 3, column = 2)

                self.logx_button = Tk.Button(
                    x_controls_frame,
                    text = 'log10',
                    command = self.toggle_log_x,
                    relief = Tk.SUNKEN if self.logx else Tk.RAISED
                )
                self.age_transform_entry = Tk.Entry(
                    x_controls_frame,
                )
                self.age_transform_entry.insert(
                    0,
                    't * (1.0 + (t / 5.0) * m**5 * 10.0**(-0.2*feh))'
                    '* m**2.3 * 10.0**(-0.4*feh)')
                self.logx_button.grid(row = 0, column = 0)
                self.age_transform_entry.grid(row = 0, column = 1)

            def create_y_controls() :
                """Create the controls for the y axis."""

                y_controls_frame = Tk.Frame(window)
                y_controls_frame.grid(row = 1, column = 1, rowspan = 2)

                self.selected_plot_quantity = Tk.StringVar()
                self.selected_plot_quantity.set(self.plot_quantities[0])
                self.plot_quantity_menu = Tk.OptionMenu(
                    y_controls_frame,
                    self.selected_plot_quantity,
                    *self.plot_quantities,
                    command = self.change_plot_quantity
                )
                self.logy_button = Tk.Button(
                    y_controls_frame,
                    text = 'log10',
                    command = self.toggle_log_y,
                    relief = Tk.SUNKEN if self.logy else Tk.RAISED
                )
                self.first_deriv_button = Tk.Button(
                    y_controls_frame,
                    text = 'd/dt',
                    command = functools.partial(self.toggle_deriv, 1)
                )
                self.second_deriv_button = Tk.Button(
                    y_controls_frame,
                    text = 'd/dt',
                    command = functools.partial(self.toggle_deriv, 2)
                )
                self.logy_button.grid(row = 0, column = 0)
                self.first_deriv_button.grid(row = 1, column = 0)
                self.second_deriv_button.grid(row = 2, column = 0)
                self.plot_quantity_menu.grid(row = 3, column = 0)

            create_x_controls()
            create_y_controls()
            Tk.Button(
                window,
                text = 'Auto Axes',
                command = self.auto_axes,
                relief = Tk.RAISED
            ).grid(row = 1,
                   column = 3,
                   sticky = Tk.N + Tk.S + Tk.W + Tk.E)

        def create_interp_controls(interp_control_frame) :
            """Create the controls for the interpolation mass & [Fe/H]."""

            self.coarse_interp_scale = dict()
            self.fine_interp_scale = dict()
            for index, quantity in enumerate(['mass', 'metallicity']) :
                self.coarse_interp_scale[quantity] = Tk.Scale(
                    interp_control_frame,
                    from_ = getattr(self, 'track_' + quantity)[0],
                    to = getattr(self, 'track_' + quantity)[-1],
                    resolution = -1,
                    length = 1000,
                    orient = Tk.HORIZONTAL,
                    command = functools.partial(self.change_interp,
                                                quantity),
                    digits = 6
                )

                self.fine_interp_scale[quantity] = Tk.Scale(
                    interp_control_frame,
                    resolution = -1,
                    length = 1000,
                    orient = Tk.HORIZONTAL,
                    command = functools.partial(self.change_interp,
                                                quantity),
                    digits = 6
                )

                self.coarse_interp_scale[quantity].grid(row = 3 * index,
                                                        column = 1)
                self.fine_interp_scale[quantity].grid(row = 3 * index + 1,
                                                      column = 1)
            Tk.Label(interp_control_frame,
                     text = 'M*/Msun').grid(row = 0, column = 0, rowspan = 2)
            Tk.Label(interp_control_frame,
                     text = '[Fe/H]').grid(row = 3, column = 0, rowspan = 2)
            ttk.Separator(interp_control_frame,
                          orient = Tk.HORIZONTAL).grid(row = 2,
                                                       column = 0,
                                                       columnspan = 2,
                                                       sticky = "ew")

        def create_track_selectors(track_selectors_frame) :
            """Create buttons to enable/disable tracks to display."""

            self.all_tracks_button = Tk.Button(
                track_selectors_frame,
                text = 'All',
                command = self.toggle_all_tracks,
                relief = Tk.RAISED
            )
            self.all_tracks_button.grid(row = 0,
                                        column = 0,
                                        rowspan = 2,
                                        columnspan = 2)

            Tk.Label(
                track_selectors_frame,
                text = 'M*/Msun'
            ).grid(row = 1,
                   column = 0,
                   columnspan = len(self.track_mass))
            Tk.Label(
                track_selectors_frame,
                text = '[Fe/H]',
            ).grid(row = 0,
                   column = 1,
                   rowspan  = len(self.track_metallicity))

            self.track_state = dict()
            self.track_mass_button = dict()
            self.track_metallicity_button = dict()
            for row, mass in enumerate(self.track_mass) :
                self.track_state[mass] = dict()
                self.track_mass_button[mass] = Tk.Button(
                    track_selectors_frame,
                    text = '%.3f' % mass,
                    command = functools.partial(self.toggle_track_mass,
                                                mass),
                    relief = Tk.RAISED
                )
                self.track_mass_button[mass].grid(row = 2 + row,
                                                  column = 1)
                for column, metallicity in enumerate(
                        self.track_metallicity
                ) :
                    if row == 0 :
                        self.track_metallicity_button[metallicity] = \
                            Tk.Button(
                                track_selectors_frame,
                                text = '%.3f' % metallicity,
                                command = functools.partial(
                                    self.toggle_track_metallicity,
                                    metallicity
                                ),
                                relief = Tk.RAISED
                            )
                        self.track_metallicity_button[metallicity].grid(
                            row = 1,
                            column = 2 + column
                        )
                    self.track_state[mass][metallicity] = Tk.IntVar()
                    Tk.Checkbutton(
                        track_selectors_frame,
                        text = '',
                        variable = self.track_state[mass][metallicity],
                        command = self.display
                    ).grid(row = 2 + row, column = 2 + column)

        def create_curve_controls(curve_control_frame) :
            """Create controls for modifying plotting."""

            Tk.Button(curve_control_frame,
                      text = 'Replot',
                      command = self.display).grid(row = 0,
                                                   column = 0,
                                                   columnspan = 2)
            curve_setup_frame = Tk.Frame(curve_control_frame)
            curve_setup_frame.grid(row = 1, column = 0)
            Tk.Label(curve_setup_frame,
                     text = 'Resolution:').grid(row = 1, column = 0)
            Tk.Entry(curve_setup_frame,
                     textvariable = self.resolution_text).grid(row = 1,
                                                               column = 1)
            Tk.Label(curve_setup_frame,
                     text = 'Track style:').grid(row = 2, column = 0)
            Tk.Entry(curve_setup_frame,
                     textvariable = self.track_plot_style).grid(row = 2,
                                                                column = 1)

            Tk.Label(curve_setup_frame,
                     text = 'Interp style:').grid(row = 3, column = 0)
            Tk.Entry(curve_setup_frame,
                     textvariable = self.interp_plot_style).grid(row = 3,
                                                                 column = 1)

            isolate_frame = Tk.Frame(curve_control_frame)
            isolate_frame.grid(row = 1, column = 1)
            Tk.Label(isolate_frame,
                     text = 'Separate Window').grid(row = 0, column = 0)
            Tk.Button(
                isolate_frame,
                text = 'Main curve',
                command = functools.partial(self.separate_window, 0)
            ).grid(row = 1, column = 0)
            Tk.Button(
                isolate_frame,
                text = 'First deriv',
                command = functools.partial(self.separate_window, 1)
            ).grid(row = 2, column = 0)
            Tk.Button(
                isolate_frame,
                text = 'Second deriv',
                command = functools.partial(self.separate_window, 2)
            ).grid(row = 3, column = 0)

        def create_main_canvas(plot_frame) :
            """Create the canvas for plotting undifferentiated quantities."""

            self.main_canvas = FigureCanvasTkAgg(self.main_figure,
                                                 master = plot_frame)
            self.main_canvas.show()
            self.main_canvas.get_tk_widget().pack(side = Tk.TOP,
                                                  fill = Tk.BOTH,
                                                  expand = 1)
            self.toolbar = NavigationToolbar2TkAgg(self.main_canvas,
                                                   plot_frame)
            self.toolbar.update()
            self.main_canvas._tkcanvas.pack(side=Tk.TOP,
                                            fill=Tk.BOTH, expand=1)
            self.main_canvas.mpl_connect('key_press_event',
                                         self.on_key_event)

        def set_initial_state() :
            """Set initial states for member variables."""

            self.main_auto_axes = True
            self.first_deriv_auto_axes = True
            self.second_deriv_auto_axes = True
            self.do_not_display = True
            self.display_job = None
            self.plot_quantities = [
                'Iconv', 'Irad', 'I', 'R', 'Lum', 'Rrad', 'Mrad'
            ]
            self.max_deriv = dict(Iconv = 2,
                                  Irad = 2,
                                  I = 2,
                                  R = 1,
                                  Lum = 0,
                                  Rrad = 1, 
                                  Mrad = 2)
            self.logx = True
            self.logy = True
            self.interpolate = False
            self.deriv = 0
            self.tracks = read_MESA(tracks_dir)
            self.track_mass = sorted(self.tracks.keys())
            self.track_metallicity = set()
            for mass_tracks in self.tracks.values() :
                self.track_metallicity.update(mass_tracks.keys())
            self.track_metallicity = sorted(list(self.track_metallicity))
            self.interp = dict(mass = 1.0, metallicity = 0.0)
            self.enabled_tracks = [
                [False for feh in self.track_metallicity]
                for m in self.track_mass
            ]

            self.resolution_text = Tk.StringVar()
            self.resolution_text.set('100')
            self.interp_plot_style = Tk.StringVar()
            self.interp_plot_style.set('.r')
            self.track_plot_style = Tk.StringVar()
            self.track_plot_style.set('xk')

            self.changing = dict(mass = False, metallicity = False)

        def configure_window() :
            """Arrange the various application elements."""

            self.window = window

            plot_frame = Tk.Frame(window)
            plot_frame.grid(row = 1,
                            column = 2,
                            rowspan = 2,
                            sticky = Tk.N + Tk.S + Tk.W + Tk.E)

            interp_control_frame = Tk.Frame(window)
            interp_control_frame.grid(row = 0, column = 1, columnspan = 3)

            plot_control_frame = Tk.Frame(window)
            plot_control_frame.grid(row = 2, column = 3)
            track_selectors_frame = Tk.Frame(plot_control_frame)
            track_selectors_frame.grid(row = 0, column = 0)
            curve_control_frame = Tk.Frame(plot_control_frame)
            curve_control_frame.grid(row = 1, column = 0)

            Tk.Grid.columnconfigure(window, 2, weight = 1)
            Tk.Grid.rowconfigure(window, 1, weight = 1)

            create_main_axes()
            create_main_canvas(plot_frame)
            create_axes_controls()
            create_interp_controls(interp_control_frame)
            create_track_selectors(track_selectors_frame)
            create_curve_controls(curve_control_frame)

            self.interpolate_button = Tk.Button(
                window,
                text = 'Interpolate',
                command = self.toggle_interpolation,
                relief = Tk.SUNKEN if self.interpolate else Tk.RAISED
            )
            self.interpolate_button.grid(row = 0,
                                         column = 0,
                                         sticky = Tk.N + Tk.S + Tk.W + Tk.E)

            interpolator_manager_frame = Tk.Frame(window)
            interpolator_manager_frame.grid(row = 1, column = 0, rowspan = 2)
            self.interpolator_manager = InterpolatorManagerGUI(
                interpolator_manager_frame,
                serialized_interpolator_dir
            )

            self.track_below = dict()

        set_initial_state()
        configure_window()

        self.change_plot_quantity(self.plot_quantities[0])
        self.change_interp('mass', self.interp['mass'])
        self.change_interp('metallicity', self.interp['metallicity'])
        self.do_not_display = False
        self.display()

def read_YREC(dirname) :
    """
    Read all YREC tracks from the given directory.
    
    Args:
        - dirname: The directory containing the tracks to read. The tracks
                   are assumed to be all files with a name ending in
                   '.track'.
    Returns: A dictionary indexed by mass (up to 3 significant figures)
             containing dictionaries indexed by metallicity (up to 3
             significant figures) each containing a dictionary with keys: t,
             Iconv, Irad, I, R, Lum, Rrad, Mrad containing the corresponding
             quantity in the units output by poet.
    """

    model_param_rex = re.compile(
        '^#Version=[ 0-9]* '
        'Mtot/Msun = +(?P<MASS>[0-9.E+-]+) +'
        'Initial: X = +[0-9.E+-]+ +Z = +[0-9.E+-]+ +'
        'Mix. length = +[0-9.E+-]+ *$'
    )
    track_fnames = glob(os.path.join(dirname, '*.track'))
    result = dict()
    Inorm = (astropy.constants.M_sun.to(astropy.units.g).value
             *
             astropy.constants.R_sun.to(astropy.units.cm).value**2)
    metallicity_key = 0.0
    all_metallicities = [metallicity_key]
    for fname in track_fnames :
        with open(fname, 'rb') as track_file :
            match = False
            while not match :
                match = model_param_rex.match(track_file.readline().decode())
            model_mass = float(match.group('MASS'))
            column_names_line = track_file.readline().decode()
            assert(column_names_line[0] == '#')
            column_names = column_names_line[1:].strip().split()
            model_data = numpy.genfromtxt(
                track_file,
                names = column_names,
                deletechars = ''
            )
            mass_key = round(model_mass, 3)
            if mass_key not in result : result[mass_key] = dict()
            assert(metallicity_key not in result[mass_key])
            result[mass_key][metallicity_key] = dict(
                t = model_data['Age(Gyr)'],
                Iconv = model_data['I_env'] / Inorm,
                Irad = model_data['I_rad'] / Inorm,
                R = 10.0**model_data['log(R/Rsun)'],
                Lum = 10.0**model_data['log(L/Lsun)'],
                Mrad = model_data['m_envp/M'] * model_mass
            )
            result[mass_key][metallicity_key]['I'] = (
                result[mass_key]['Iconv']
                +
                result[mass_key]['Irad']
            )
            Rrad_key = ('r_envp/M' if 'r_envp/M' in model_data.dtype.names
                        else 'r_envp/R')
            result[mass_key][metallicity_key]['Rrad'] = (
                model_data[Rrad_key]
                * 
                result[mass_key]['R']
            )
    return result

def read_MESA(dirname) :
    """
    Read all MESA tracks from the given directory.
    
    Args:
        - dirname: The directory containing the tracks to read. The tracks
                   are assumed to be all files with filenames like
                   M<stellar_mass>_Z<metallicity>.csv
    Returns: A dictionary indexed by mass (up to 3 significant figures)
             containing dictionaries indexed by metallicity (up to 3
             significant figures) each containing a numpy record array with
             keys: t, Iconv, Irad, I, R, Lum, Rrad, Mrad containing the
             corresponding quantity in the units output by poet.
    """

    def translate_name(name) :
        """Return the output record name for the given name in the track."""

        if name == 'age' : return 't'
        elif name == 'R_star' : return 'R'
        elif name == 'L_star' : return 'Lum'
        elif name == 'R_tachocline' : return 'Rrad'
        else : return name.replace('_', '')

    result = dict()
    fname_rex = re.compile(
        'M(?P<MASS>[0-9.E+-]+)_Z(?P<METALLICITY>[0-9.E+-]+).csv'
    )
    track_fnames = glob(os.path.join(dirname, '*.csv'))
    for fname in track_fnames :
        parsed_fname = fname_rex.match(os.path.basename(fname))
        if not parsed_fname : 
            print('Skipping ' + repr(fname))
            continue
        mass_key = round(float(parsed_fname.group('MASS')), 3)
        feh = stellar_evolution_library.feh_from_z(
            float(parsed_fname.group('METALLICITY'))
        )
        metallicity_key = round(
            (-1.0 if feh < 0 else 1.0) * floor(abs(feh) * 1000) / 1000,
            3
        )
        if mass_key not in result : result[mass_key] = dict()
        assert(metallicity_key not in result[mass_key])
        result[mass_key][metallicity_key] = numpy.genfromtxt(fname,
                                                             delimiter = ',',
                                                             names = True)
        result[mass_key][metallicity_key].dtype.names = tuple(
            translate_name(name) for name in
            result[mass_key][metallicity_key].dtype.names
        )
        result[mass_key][metallicity_key]['t'] /= 1e9
    return result

if __name__ == '__main__' :
    main_window = Tk.Tk()
    main_window.wm_title("Stellar Evolution Interpolation Explorer")
    ap = InterpolationInteractive(main_window, mesa_track_dir)

    Tk.mainloop()
