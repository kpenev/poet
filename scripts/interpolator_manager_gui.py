#!/usr/bin/python3

import sys
sys.path.append('../PythonPackage')
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.library_interface import MESAInterpolator
from GUIUtil import Dialog

if sys.version_info[0] < 3:
    import Tkinter as Tk
    from Tkinter import ttk
else:
    import tkinter as Tk
    from tkinter import ttk

class NewInterpolatorDialog(Dialog) :
    """A dialogue that popus up for new interpolator generation."""

    def __init__(self, parent, forbidden_names) :
        """Set up the dialog."""

        super().__init__(parent = parent, title = 'Create New Interpolator')
        self.forbidden_names = forbidden_names
        self.create = False

    def body(self, master) :
        """Add entry for name and OK and Cancel buttons."""

        Tk.Label(master, text = 'Create new iterator').grid(row = 0,
                                                            column = 0,
                                                            columnspan = 2)
        Tk.Label(master, text = 'Name: ').grid(row = 1, column = 0)
        self.name = Tk.StringVar()

        Tk.Entry(master, textvariable = self.name).grid(row = 1,
                                                      column = 1,
                                                      sticky = Tk.E + Tk.W)

    def apply(self) :
        """Note that the new interpolator should be generated."""

        self.create = True

    def validate(self) :
        return (self.name.get()
                and
                self.name.get() not in self.forbidden_names)

class InterpolatorManagerGUI(StellarEvolutionManager) :

    def _add_selector(self, parent, choices, tk_variable, command) :
        """Add a drop-down menu for the list of available interpolators."""

        menu = Tk.OptionMenu(parent,
                             tk_variable,
                             *choices,
                             command = command)
        menu.grid(row = 0, column = 0)
        menu.grid_anchor(Tk.N)

    def _add_track_selectors(self, orientation = 'vertical') :
        """Create check buttons indicating tracks an interpolator uses."""

        for widget in self.track_selector_parent.winfo_children() :
            widget.destroy()

        tracks = self.get_suite_tracks(self.selected_suite.get())
        masses = sorted({t.mass for t in tracks})
        metallicities = sorted({t.metallicity for t in tracks})
        self.mass_selected = {m: Tk.IntVar() for m in masses}
        self.metallicity_selected = {feh: Tk.IntVar()
                                     for feh in metallicities}

        row, column = 1, 1
        Tk.Label(self.track_selector_parent, text = 'M*').grid(
            row = row, column = column
        )
        for m in masses :
            if orientation == 'horizontal' : column += 1
            else : row += 1
            Tk.Checkbutton(
                self.track_selector_parent,
                text = str(m),
                variable = self.mass_selected[m]
            ).grid(row = row, column = column)

        if orientation == 'horizontal' : row, column = 2, 1
        else : column, row = 2, 1
        Tk.Label(self.track_selector_parent, text = '[Fe/H]').grid(
            row = row, column = column
        )
        for feh in metallicities :
            if orientation == 'horizontal' : column += 1
            else : row += 1
            Tk.Checkbutton(
                self.track_selector_parent,
                text = str(feh),
                variable = self.metallicity_selected[feh]
            ).grid(row = row, column = column)

    def _add_interp_param_controls(self, parent, orientation = 'vertical') :
        """Create the controls for displaying the selected interpolator."""

        self.nodes = {
            quantity: Tk.StringVar()
            for quantity in MESAInterpolator.quantity_list
        }
        self.smoothing = {
            quantity: Tk.StringVar()
            for quantity in MESAInterpolator.quantity_list
        }

        for quantity_name, quantity_index in \
                MESAInterpolator.quantity_ids.items() :
            if orientation == 'horizontal' :
                row, column = 1, quantity_index + 2
            else :
                assert(orientation == 'vertical')
                column, row = 1, quantity_index + 2
            Tk.Label(parent,
                     text = quantity_name).grid(row = row,
                                                column = column)
            if orientation == 'horizontal' : row += 1
            else : column += 1
            Tk.Entry(
                parent,
                textvariable = self.nodes[quantity_name],
                width = 10
            ).grid(row = row, column = column)
            if orientation == 'horizontal' : row += 1
            else : column += 1
            Tk.Entry(
                parent,
                textvariable = self.smoothing[quantity_name],
                width = 10
            ).grid(row = row, column = column)

        if orientation == 'horizontal' : row, column = 2, 1
        else : column, row = 2, 1
        Tk.Label(parent,
                 text = 'Nodes:').grid(row = row, column = column)
        if orientation == 'horizontal' : row += 1
        else : column += 1
        Tk.Label(parent,
                 text = 'Smoothing:').grid(row = row, column = column)

    def _refresh_interpolator(self, *ignore) :
        """
        Update the display after an interpolator change.

        Args:
            Ignored, but allowed to accomodate tkinter generated calls.

        Returns: None
        """
        
        self.interpolator = self.get_interpolator_by_name(
            self.selected_interpolator.get()
        )
        for quantity in self.interpolator.quantity_list :
            self.nodes[quantity].set(str(self.interpolator.nodes[quantity]))
            self.smoothing[quantity].set(
                str(self.interpolator.smoothing[quantity])
            )
        self.selected_suite.set(self.interpolator.suite)
        self._refresh_suite()
        for mass in self.interpolator.track_masses :
            self.mass_selected[mass].set(1)
        for metallicity in self.interpolator.track_metallicities :
            self.metallicity_selected[metallicity].set(1)

    def _refresh_suite(self, *ignore) :
        """
        Update the display after a suite change.

        Args:
            Ignored, but allowed to accomodate tkinter generated calls.

        Returns: None
        """
        
        self._add_track_selectors(orientation = self.orientation)
        for tkvar in self.mass_selected.values() : tkvar.set(0)
        for tkvar in self.metallicity_selected.values() : tkvar.set(0)

    def _match_config(self) :
        """
        Set the current interpolator per current setup (create if necessary).
        """

        masses = []
        for mass, enabled in sorted(self.mass_selected.items()) :
            if enabled.get() : masses.append(mass)

        metallicities = []
        for metallicity, enabled in sorted(
                self.metallicity_selected.items()
        ) :
            if enabled.get() : metallicities.append(metallicity)

        nodes = {quantity: int(tkvar.get())
                 for quantity, tkvar in self.nodes.items()}

        smoothing = {quantity: float(tkvar.get())
                     for quantity, tkvar in self.smoothing.items()}

        self.interpolator = self.get_interpolator(
            masses = masses,
            metallicities = metallicities,
            model_suite = self.selected_suite.get(),
            nodes = nodes,
            smoothing = smoothing
        )
        if self.interpolator is None :
            dialog = NewInterpolatorDialog(self.parent, self._interp_list)
            self.parent.wait_window(dialog)
            if dialog.create :
                self.interpolator = self.get_interpolator(
                    masses = masses,
                    metallicities = metallicities,
                    model_suite = self.selected_suite.get(),
                    nodes = nodes,
                    smoothing = smoothing,
                    new_interp_name = dialog.name.get()
                )
                self._interp_list.append(self.interpolator.name)
                self._add_selector(self._interp_selector_frame,
                                   self._interp_list,
                                   self.selected_interpolator,
                                   self._refresh_interpolator)
        self.selected_interpolator.set(self.interpolator.name)
        self._refresh_interpolator()

    def __init__(self,
                 parent,
                 serialization_path,
                 orientation = 'vertical') :
        """
        Create a manager and set-up its GUI in a parent.

        Args:
            - parent:
                A Tk object under which to place all widgets of the manager.
            - serialization_path:
                The path where to store serialized interpolators.
            - orientation:
                Should the widget arrangement span more horizontally or
                vertically. Must be have a value of either 'horizontal' or
                'vertical'.

        Returns: None
        """

        super().__init__(serialization_path)

        self.orientation = orientation
        self.parent = parent

        self._interp_list = self.list_interpolator_names()
        self.selected_interpolator = Tk.StringVar()
        self.selected_interpolator.set(self._interp_list[0])
        self._interp_selector_frame = Tk.Frame(parent)
        self._interp_selector_frame.grid_anchor(Tk.N)
        self._interp_selector_frame.grid(row = 0, column = 0)
        self._add_selector(self._interp_selector_frame,
                           self._interp_list,
                           self.selected_interpolator,
                           self._refresh_interpolator)

        interp_param_frame = Tk.Frame(parent)
        interp_param_frame.grid(row = 1, column = 0, sticky = Tk.W + Tk.E)
        self._add_interp_param_controls(interp_param_frame,
                                        orientation = orientation)

        suite_list = self.list_suites()
        self.selected_suite = Tk.StringVar()
        self.selected_suite.set(suite_list[0])
        suite_selector_frame = Tk.Frame(parent)
        suite_selector_frame.grid_anchor(Tk.CENTER)
        suite_selector_frame.grid(row = 2, column = 0)
        self._add_selector(suite_selector_frame,
                           suite_list,
                           self.selected_suite,
                           self._refresh_suite)

        self.track_selector_parent = Tk.Frame(parent)
        self.track_selector_parent.grid(row = 3, column = 0)
        get_interp_frame = Tk.Frame(parent)
        get_interp_frame.grid(row = 4, column = 0)
        Tk.Button(get_interp_frame,
                  text = 'Get Interpolator',
                  command = self._match_config).grid(row = 0, column = 0)
        self._refresh_interpolator()

    def current_interpolator(self) :
        """Return the currently selected interpolator."""

        return self.interpolator

if __name__ == '__main__' :
    main_window = Tk.Tk()
    main_window.wm_title("Interpolator Manager")
    ap = InterpolatorManagerGUI(main_window,
                                '../stellar_evolution_interpolators')

    Tk.mainloop()
