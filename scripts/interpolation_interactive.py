#!/usr/bin/python3

import matplotlib

matplotlib.use('TkAgg')

import sys

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
import subprocess
import psutil
import functools
import time

if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

class Structure :
    """An empty class used only to hold user defined attributes."""

    def __init__(self, **initial_attributes) :
        """Create a class with (optionally) initial attributes."""

        for attribute_name, attribute_value in initial_attributes.items() :
            setattr(self, attribute_name, attribute_value)

    def format(self, prefix='') :
        """Generate a tree-like representation of self."""

        result=''
        for attr_name in dir(self) :
            if attr_name[0]!='_' :
                attribute=getattr(self, attr_name)
                if isinstance(attribute, Structure) :
                    result+=(prefix
                             +
                             '|-'
                             +
                             attr_name
                             +
                             '\n'
                             +
                             attribute.format(prefix + '| '))
                else : result+=(prefix
                                +
                                '|-'
                                +
                                attr_name
                                +
                                ': '
                                +
                                str(attribute)
                                +
                                '\n')
        return result

def star_lifetime(mass) : return min(9.0 * mass**-3, 10.0)

class PoetInterp :
    """Interface with POET to perform stellar evolution interpolation."""

    def start_poet(self) :
        """(Re-)start the poet process used for interpolation."""

        self.poet = subprocess.Popen(
            [
                self.poet_executable,
                '--input-columns', ','.join(self.input_columns),
                '--output-columns', ','.join(self.output_columns)
            ],
            stdin = subprocess.PIPE,
            bufsize = 0
        )
        self.poet_process_files = psutil.Process(self.poet.pid).open_files

    def __init__(self, poet_executable, output_columns) :
        """Start poet and wait for interpolation requests."""

        self.input_columns = ['M', 't0', 'tmax', 'maxdt', 'outf']
        self.input_line_format = ' '.join(
            ['%(' + c + ')s' for c in self.input_columns]
        ) + '\n'
        self.output_columns = output_columns
        self.poet_executable = poet_executable
        self.start_poet()

        self.output_fname = 'poet.evol'

    def __call__(self, star_mass, start_age, end_age, max_age_step) :
        """Return the interpolated stellar evolution per the arguments."""

        assert(not os.path.exists(self.output_fname))
        end_age = min(end_age, star_lifetime(star_mass))
        input_line = (self.input_line_format
                      %
                      dict(M = repr(star_mass),
                           t0 = repr(start_age),
                           tmax = repr(end_age),
                           maxdt = repr(max_age_step),
                           outf = self.output_fname)).encode('ascii')

        print(input_line)
        self.poet.stdin.write(input_line)
        while(not os.path.exists(self.output_fname)) : time.sleep(0.01)

        last_age = numpy.nan
        while not (end_age - last_age < 0.5 * max_age_step) :
            if self.poet.poll() is not None :
                print('Restarting poet!')
                self.start_poet()
                break
            try :
                result = numpy.genfromtxt(self.output_fname,
                                          names = True)
                last_age = result['t'][-1]
            except : 
                last_age = numpy.nan
                if self.poet.poll() is not None :
                    print('Restarting poet!')
                    self.start_poet()
                    self.poet.stdin.write(input_line)

        os.remove(self.output_fname)
        return result

class Application :

    def quit(self) :
        """Exit the application."""

        main_window.quit()     # stops mainloop
        main_window.destroy()  # this is necessary on Windows to prevent
                               # Fatal Python Error: PyEval_RestoreThread:
                               # NULL tstate

    def display(self) :
        """(Re-)draw the plot as currently configured by the user."""

        def plot_interpolation(star_mass, plot_func) :
            """Plot an interpolated stellar evolution track."""

            interpolated = self.poet_interp(star_mass,
                                            0,
                                            star_lifetime(star_mass),
                                            0.001)
            plot_func(interpolated['t'],
                      interpolated[self.plot_quantity],
                      '.r')

            if self.deriv > 0 :
                d1_y = numpy.copy(interpolated['D' + self.plot_quantity])
                if self.deriv > 1 :
                    d2_y = numpy.copy(interpolated['DD' + self.plot_quantity])

                if self.logy : 
                    d1_y /= interpolated[self.plot_quantity]
                    if self.deriv > 1 :
                        d2_y = (d2_y / interpolated[self.plot_quantity]
                                -
                                d1_y**2)

                if self.logx : 
                    deriv_plot_funcname = 'semilogx'
                    d1_y *= interpolated['t']
                    if self.deriv > 1 :
                        d2_y = d1_y + interpolated['t']**2 * d2_y
                else :
                    deriv_plot_funcname = 'plot'


                getattr(self.first_deriv_axes, deriv_plot_funcname)(
                    interpolated['t'],
                    d1_y,
                    '.r'
                )
                if self.deriv > 1 :
                    getattr(self.second_deriv_axes, deriv_plot_funcname)(
                        interpolated['t'],
                        d2_y,
                        '.r'
                    )
        if self.do_not_display : return
        
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

        plot_interpolation(self.interp_mass, plot)

        for track_index in range(len(self.track_masses)) :
            if (
                    self.enabled_tracks[track_index]
                    or
                    track_index == self.track_index_below
                    or
                    track_index == self.track_index_below + 1
            ) :
                track_mass = self.track_masses[track_index]

                if self.enabled_tracks[track_index] :
                    plot_interpolation(track_mass, plot)

                plot(self.tracks[track_mass]['t'],
                     self.tracks[track_mass][self.plot_quantity],
                     'xk')
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
                    self.second_deriv_axes.set_xlim(self.main_axes.get_xlim())
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
        self.display()

    def change_interp_mass(self, new_mass) :
        """Modify the stellar mass for which to display the interpolation."""

        if self.changing_mass : return
        self.changing_mass = True

        self.interp_mass = float(new_mass)
        self.track_index_below = len(self.track_masses) - 2
        while (
                self.track_index_below >= 0
                and
                self.track_masses[self.track_index_below] > \
                self.interp_mass
        ) : self.track_index_below -= 1

        self.fine_interp_mass_scale.config(
            from_ = (self.track_masses[self.track_index_below]
                     +
                     1e-5),
            to = (self.track_masses[self.track_index_below + 1]
                  -
                  1e-5)
        )

        self.coarse_interp_mass_scale.set(self.interp_mass)
        self.fine_interp_mass_scale.set(self.interp_mass)

        if self.display_job : self.main_window.after_cancel(self.display_job)
        self.display_job = self.main_window.after(100, self.display)

        self.changing_mass = False
        return True

    def toggle_track(self, track_index) :
        """Toggle displaying the track with the given index."""

        self.enabled_tracks[track_index] = (
            not self.enabled_tracks[track_index]
        )
        self.track_buttons[track_index].config(
            relief = (Tk.SUNKEN if self.enabled_tracks[track_index]
                      else Tk.RAISED)
        )
        self.display()

    def set_all_tracks(self, state) :
        """Enable displaying all tracks."""

        for i in range(len(self.enabled_tracks)) :
            self.enabled_tracks[i] = state
            self.track_buttons[i].config(
                relief = Tk.SUNKEN if state else Tk.RAISED
            )

        self.display()

    def auto_axes(self) :
        """Re-plot letting matplotlib determine the axes limits."""

        self.main_auto_axes = True
        self.first_deriv_auto_axes = True
        self.second_deriv_auto_axes = True
        self.display()

    def __init__(self, main_window, tracks) :
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

            y_controls_frame = Tk.Frame(main_window)
            y_controls_frame.grid(row = 1, column = 0)

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

            self.logx_button = Tk.Button(
                main_window,
                text = 'log10',
                command = self.toggle_log_x,
                relief = Tk.SUNKEN if self.logx else Tk.RAISED
            )
            self.logx_button.grid(row = 2, column = 1)

            Tk.Button(
                main_window,
                text = 'Auto Axes',
                command = self.auto_axes,
                relief = Tk.RAISED
            ).grid(row = 0,
                   column = 2,
                   sticky = Tk.N + Tk.S + Tk.W + Tk.E)

        def create_mass_controls(mass_control_frame) :
            """Create the controls to select the interpolation mass."""

            self.coarse_interp_mass_scale = Tk.Scale(
                mass_control_frame,
                from_ = self.track_masses[0],
                to = self.track_masses[-1],
                resolution = -1,
                length = 1000,
                orient = Tk.HORIZONTAL,
                command = self.change_interp_mass,
                digits = 6
            )

            self.fine_interp_mass_scale = Tk.Scale(
                mass_control_frame,
                resolution = -1,
                length = 1000,
                orient = Tk.HORIZONTAL,
                command = self.change_interp_mass,
                digits = 6
            )

            self.coarse_interp_mass_scale.grid(row = 0, column = 0)
            self.fine_interp_mass_scale.grid(row = 1, column = 0)

        def create_track_selectors(track_selectors_frame) :
            """Create buttons to enable/disable tracks to display."""

            self.track_buttons = [
                Tk.Button(
                    track_selectors_frame,
                    text = "M=%.3f" % mass,
                    command = functools.partial(self.toggle_track, index)
                )
                for index, mass in enumerate(self.track_masses)
            ]

            Tk.Button(
                track_selectors_frame,
                text = 'All',
                command = functools.partial(self.set_all_tracks, True)
            ).grid(row = 0, column = 0)
            Tk.Button(
                track_selectors_frame,
                text = 'None',
                command = functools.partial(self.set_all_tracks, False)
            ).grid(row = 0, column = 1)

            for index, button in enumerate(self.track_buttons) :
                button.grid(row = index + 1, column = 0, columnspan = 2)

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
                              Rrad = 2, 
                              Mrad = 1)
        self.logx = False
        self.logy = False
        self.deriv = 0
        self.tracks = tracks
        self.track_masses = sorted(tracks.keys())
        self.interp_mass = 1.0
        self.enabled_tracks = [False for m in self.track_masses]
        self.changing_mass = False
        self.main_window = main_window

        plot_frame = Tk.Frame(main_window)
        plot_frame.grid(row = 1,
                        column = 1,
                        sticky = Tk.N + Tk.S + Tk.W + Tk.E )

        mass_control_frame = Tk.Frame(main_window)
        mass_control_frame.grid(row = 0, column = 1)

        track_selectors_frame = Tk.Frame(main_window)
        track_selectors_frame.grid(row = 1, column = 2)

        Tk.Grid.columnconfigure(main_window, 1, weight = 1)
        Tk.Grid.rowconfigure(main_window, 1, weight = 1)

        create_main_axes()
        create_main_canvas(plot_frame)
        create_axes_controls()
        create_mass_controls(mass_control_frame)
        create_track_selectors(track_selectors_frame)

        self.change_interp_mass(self.interp_mass)
        self.change_plot_quantity(self.plot_quantities[0])
        self.do_not_display = False
        poet_output_columns = ['t']
        for quantity, max_deriv in self.max_deriv.items() :
            poet_output_columns.append(quantity)
            if max_deriv > 0 : poet_output_columns.append('D' + quantity)
            if max_deriv > 1 : poet_output_columns.append('DD' + quantity)
        self.poet_interp = PoetInterp('./poet', poet_output_columns)
        self.display()

def read_YREC(dirname) :
    """
    Read all YREC tracks from the given directory.
    
    Args:
        - dirname: The directory containing the tracks to read. The tracks
                   are assumed to be all files with a name ending in
                   '.track'.
    Returns: A dictionary with keys the track masses (to 3 significant
             figures), each containing a dictionary with keys:
             t, Iconv, Irad, I, R, Lum, Rrad, Mrad containing the
             corresponding quantity in the units output by poet.
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
            result[mass_key] = dict(
                t = model_data['Age(Gyr)'],
                Iconv = model_data['I_env'] / Inorm,
                Irad = model_data['I_rad'] / Inorm,
                R = 10.0**model_data['log(R/Rsun)'],
                Lum = 10.0**model_data['log(L/Lsun)'],
                Mrad = model_data['m_envp/M'] * model_mass
            )
            result[mass_key]['I'] = (result[mass_key]['Iconv']
                                     +
                                     result[mass_key]['Irad'])
            Rrad_key = ('r_envp/M' if 'r_envp/M' in model_data.dtype.names
                        else 'r_envp/R')
            result[mass_key]['Rrad'] = (model_data[Rrad_key]
                                        * 
                                        result[mass_key]['R'])
    return result

if __name__ == '__main__' :
    tracks = read_YREC('../poet_src/YREC')
    main_window = Tk.Tk()
    main_window.wm_title("Stellar Evolution Interpolation Explorer")
    ap = Application(main_window, tracks)

    Tk.mainloop()
