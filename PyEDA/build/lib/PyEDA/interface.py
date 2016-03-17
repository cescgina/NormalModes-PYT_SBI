"""
Interface Module:
This module is part of the PyEDA program. It provides several classes to use
PyEDA with a graphical Interface using Tkinter. It also acts as a main script
to drive the EDA analyses.

Classes:
    EDA_app -> It is the main window where all the other frames are projected.
    inherited from tkinter.Tk.
    StartPage -> Is the first frame the users sees after the execution of PyEDA.
    It leads to the loading data frame or the "about us" page.
    inherited from tkinter.Frame
    initial_root -> Acts as a loading data frame. Is where the user interacts
    with the program by entering the required data and options.
    inherited from tkinter.Frame
    waiting_window -> Is where the analysis is taking place. Incorporates the
    same steps that the PyEDA.py CLI. Leads to the plotting frame.
    inherited from tkinter.Frame
    plot_window -> It shows buttons to the different ploting capabilities of
    the software. It integrates matplotlib in the Tkinter window.
"""

import tkinter
from tkinter import ttk
import os.path
import os
import sys
from tkinter import filedialog
import PyEDA.edanalysis as eda
import PyEDA.helper_module as mdl
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import Bio.PDB as pdb


# implement the default mpl key bindings


TITLE_FONT = ("Helvetica", 18, "bold")
text_fornt = ("Courier", 14)
###############################################################################
class EDA_app(tkinter.Tk):
    """
    It is the Tk window. It has the title information and it keeps the app.data dictionary where the
    diferent frames can share data. It also have a method to switch between
    frames. This is useful to maintain a program-like expirience.

    Methods:
        show_frame -> It's the only method in this class and it is used to
        switch to a specific frame whereever frame the user is.

    Attributes:
        self.app_data -> It is the data Hub for the whole application. Is where
        the data is kept when switching from frame to frame
        self.title -> It sets the title of the application.
        self.frames -> It is a dictionary that contains the frame instances all
        preloaded from the begining. It is used by show_frame method to
        get redirection to specific frames.
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        tkinter.Tk.__init__(self, *args, **kwargs)
        self.app_data = {"filename": '',
                         "pdbid": '',
                         "plot": '',
                         "pdbfilename": '',
                         "pathname": '',
                         "atom": '',
                         "mode": '',
                         "pathplots": '',
                         "RMSD_plot": ''}
        # the container is where we'll stack a bunch of frames
        # on top of each other, then the one we want visible
        # will be raised above the others
        self.title("PyEDA: by JF Gilabert & D Mas")
        #self.geometry("600x600")
        container = tkinter.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (StartPage, initial_root, waiting_window, About_EDA, \
            plot_window):
            page_name = F.__name__
            frame = F(container, self)
            self.frames[page_name] = frame
            frame.grid(row=0, column=0, sticky="nsew")


        self.show_frame("StartPage")

    def show_frame(self, page_name):
        """Show a frame for the given page name"""
        frame = self.frames[page_name]
        frame.tkraise()


###############################################################################
class StartPage(tkinter.Frame):
    """
    It is the Starting page. It has few information displayed about the program.
    it let the user go to the analysis or the "about us" page.

    Methods:
        No methods implemented

    Attributes:
        self.controller -> It is a reference to the main Tk window.
    """

    def __init__(self, parent, controller):
        """
        Frame Constructor. We use the controller atribute to refer to the
        Tk window that is "controling" this frame. This way it can acces
        to the main window methods.
        """
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        title = tkinter.Label(self, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25))
        title.grid(row=1, column=0, columnspan=3)
        photo = tkinter.PhotoImage(file="image.gif")
        w = tkinter.Label(self, image=photo)
        w.photo = photo
        w.grid(row=0, column=0, columnspan=3)
        subtitle = tkinter.Label(self, text="""
        This Software is the result for a project of Structual Bioinformatics
        and Introduction to Python from the Master of Science in Bioinformatics
        at Universitat Pompeu Fabra. We use the Essential Dynamics approach
        in order to analyse the Normal Modes of NMR structures and Molecular
        Dynamics. The input format is pdb although it  can be also retrieved
        automatically by the software.
        """, font=text_fornt, justify=tkinter.LEFT)
        subtitle2 = tkinter.Label(self, text="Authors: JF Gilabert & D Mas",\
         font=("Helvetica", 22, "bold"))
        subtitle2.grid(row=2, column=0, columnspan=3)
        subtitle.grid(row=3, column=0, columnspan=3)


        button1 = tkinter.Button(self, text="Start the Analysis",
                            command=lambda: controller.show_frame("initial_root"))
        button2 = tkinter.Button(self, text="About PyEDA",
                            command=lambda: controller.show_frame("About_EDA"))
        button3 = tkinter.Button(self, text="Close App", command=self.quit)
        button1.grid(row=5, column=0)
        button2.grid(row=5, column=1)
        button3.grid(row=5, column=2)



################################################################################

class initial_root(tkinter.Frame):
    """
    This is the initial page of the application where the files are loaded and
    where the options can be selected.

    Methods:
        select_a -> It shows a message with the atom chosen
        select_m -> It shows a message with the mode chosen
        check_uniprot_accession_code -> Checks if the code have len==4
        get_pdb -> pass on the pdb code and move the user to the waiting_window
        select_file -> pass on the pdb filename and move the user
        to the waiting_window

    Attributes:
        self.controller -> It is a reference to the main Tk window.
        self.entry_var -> to let the user place a pdb code
        self.labels -> to display the method chosen in the radio buttons
        self.m -> to pass on the mode chosen
        self.var -> to pass on the atoms chosen
    """
    def __init__(self, parent, controller):
        """
        Frame Constructor. We use the controller atribute to refer to the
        Tk window that is "controling" this frame. This way it can acces
        to the main window methods.
        There are diverse types of buttons used to load the data.
        """
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        ### title
        title = tkinter.Label(self, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25), justify= "center")
        title.grid(row=1, column=0, columnspan=5)
        ### photo
        photo = tkinter.PhotoImage(file="image.gif")
        w = tkinter.Label(self, image=photo)
        w.photo = photo
        w.grid(row=0, column=0, columnspan=5)
        ### subtitle
        subtitle = tkinter.Label(self, text="Enter a file or a pdb code:", \
                                 font=("Helvetica", 20))
        subtitle2 = tkinter.Label(self, text="Parameters:", \
                                 font=("Helvetica", 20))
        subtitle.grid(row=2, column=0, columnspan=2)
        subtitle2.grid(row=2, column=2, columnspan=3)
        ### pdb code entry
        self.entry_var = tkinter.StringVar()
            ### entry
        label = tkinter.Label(self, text="pdb code:",
                              font=("Helvetica", 15))
        label.grid(row=3, column=0)
        self.entry = tkinter.Entry(self, bd=2, width=5, \
            textvariable=self.entry_var)
        self.entry.grid(row=4, column=0)
            ### check button
        b = tkinter.Button(self, text="Check",
                           command=self.check_uniprot_accession_code)
        b.grid(row=5, column=0)
            ### get pdb button
        get = tkinter.Button(self, text="Get PDB", command=self.get_pdb)
        get.grid(row=6, column=0, columnspan=2)
            ### add a file button
        add = tkinter.Button(self, text="Add a pdb file", \
            command=self.select_file)
        add.grid(row=4, column=1, rowspan=2)

        label4 = tkinter.Label(self, text="Go Back:",
                              font=("Helvetica", 15))
        label4.grid(row=7, column=0, columnspan=2)
        ### choose atoms
        label2 = tkinter.Label(self, text="Mode:",
                              font=("Helvetica", 15))
        label2.grid(row=3, column=2, columnspan=3)
        self.m = tkinter.StringVar()
        m1 = tkinter.Radiobutton( self, text="NMR", variable = self.m, \
            value = 'NMR', command = self.select_m)
        m2 = tkinter.Radiobutton( self, text="MD", variable = self.m, \
            value = 'MD', command = self.select_m )

        m1.grid( row=4, column=2)
        m2.grid( row=4, column=3, columnspan=2)

        self.label_m = tkinter.Label(self)
        self.label_m.grid(row=5, column=2, columnspan=3)

        ### choose mode
        label3 = tkinter.Label(self, text="Atoms:",
                              font=("Helvetica", 15))
        label3.grid(row=6, column=2, columnspan=3)
        self.var = tkinter.StringVar()
        r1 = tkinter.Radiobutton( self, text="Only CA", variable = self.var,\
            value = 'CA', command = self.select_a)
        r2 = tkinter.Radiobutton( self, text="Backbone", \
            variable = self.var, value = 'Back', command = self.select_a )
        r3 = tkinter.Radiobutton( self, text="All", variable = self.var, \
            value = 'all', command = self.select_a )

        r1.grid( row=7, column=2)
        r2.grid( row=7, column=3)
        r3.grid( row=7, column=4)

        self.label_a = tkinter.Label(self)
        self.label_a.grid(row=8, column=2, columnspan=3)

        ### start page button
        button = tkinter.Button(self, text="Start page",
                           command=lambda: controller.show_frame("StartPage"))
        button.grid( row=8, column=0)
            ### close button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.grid( row=8, column=1)
        ### Weighting loop
        i=0
        while (i<4):
            self.grid_columnconfigure(i, weight=1)
            i+=1
    def select_m(self):
        """displays a message with the option chosen"""
        selection = "You selected the option %s" %self.m.get()
        self.label_m.config(text=selection)
        self.controller.app_data["mode"] = self.m.get()

    def select_a(self):
        """displays a message with the option chosen"""
        selection = "You selected the option %s" %self.var.get()
        self.label_a.config(text=selection)
        self.controller.app_data["atom"] = self.var.get()

    def check_uniprot_accession_code(self):
        """ uses a function to check weather the code in 4 chars long"""
        if mdl.pdb_code_check(self.entry_var.get()):
            self.entry["foreground"] = "green"
        else:
            sys.stderr.write("%s is not a uniprotaccession code format\n" \
                            %self.entry_var.get())
            self.entry["foreground"] = "red"
    def get_pdb(self):
        """gets the pdb code and stores the info, then move on to waiting"""
        interface_code = str(self.entry_var.get())
        if mdl.pdb_code_check(interface_code):
            self.controller.app_data["pdbid"] = interface_code
            self.controller.app_data["pdbfilename"] = 'pdb'+interface_code+'.ent'
            self.controller.app_data["pathname"] = 'pdbfiles/'
            self.controller.app_data["pathplots"] = 'plots/'
            if not os.path.exists(self.controller.app_data["pathname"]):
                os.mkdir(self.controller.app_data["pathname"])
            if not os.path.exists(self.controller.app_data["pathplots"]):
                os.mkdir(self.controller.app_data["pathplots"])
            self.controller.show_frame("waiting_window")

        else:
            sys.stderr.write("%s is not a uniprotaccession code format\n" \
                            %self.entry_var.get())
            self.entry["foreground"] = "red"

    def select_file(self):
        """gets the pdb file name and stores the info, then move on to waiting"""
        filename = filedialog.askopenfilename()
        if not (filename.endswith('pdb') or filename.endswith('ent')):
            raise ValueError('Your input file is not a valid PDB file, \
please use a pdb or ent file')
        path_name,pdb_file = os.path.split(filename)
        sys.stderr.write("Reading the structure from {}.\n".format(filename))
        self.controller.app_data["pathname"] = path_name+'/'
        self.controller.app_data["pathplots"] = path_name+'/plots/'
        self.controller.app_data["pdbfilename"] = pdb_file
        if pdb_file.endswith('pdb'):
            pdb_id = pdb_file[:4]
        elif pdb_file.endswith('ent'):
            pdb_id = pdb_file[3:7]
        else:
            raise ValueError('Your input file is not a valid PDB file, \
please use a pdb or ent file')
        self.controller.app_data["pdbid"] = pdb_id
        self.controller.app_data["atom"] = self.var
        self.controller.show_frame("waiting_window")


################################################################################

class About_EDA(tkinter.Frame):
    """
    It is a information page that contains an explanation of the software
    features.

    Methods:

    Attributes:

    """
    def __init__(self, parent, controller):
        """ Constructor """
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        label1 = tkinter.Label(self, text="Help and Documentation", font=TITLE_FONT)
        label2 = tkinter.Label(self, text="""
        PyEDA is a python based software that performs a Normal Mode Analysis
        using an Essential Dynamics approach.
        """, font=text_fornt, justify=tkinter.LEFT)
        label3 = tkinter.Label(self, text="""
        You can perform the analysis from a pdb file or you can introduce the
        desired pdb code and the program will retrieve the pdb file for you if
        you have internet connection.
        """, font=text_fornt, justify=tkinter.LEFT)
        label4 = tkinter.Label(self, text="""
        The Program will output a set of files and plots. First, it will
        generate a superimposed pdb files with all the NMR models or MD frames.
        After that it will calculate the covariance matrix of the coordinates
        and will output a plot of the eigenvalues derived from the
        diagonalitzation of the matrix. Finally, it will output a pdb file
        with the moved coordinates.
        """, font=text_fornt, justify=tkinter.LEFT)
        label1.grid(row=0,column=0,columnspan=3)
        label2.grid(row=1,column=0,columnspan=3)
        label3.grid(row=2,column=0,columnspan=3)
        label4.grid(row=3,column=0,columnspan=3)

        button = tkinter.Button(self, text="Go to the start page",
                           command=lambda: controller.show_frame("StartPage"))
        button.grid(row=4,column=0)

        button1 = tkinter.Button(self, text="Start the Analysis",
                            command=lambda: controller.show_frame("initial_root"))
        button1.grid(row=4,column=1)

        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.grid(row=4,column=2)
        ### Weighting loop
        i=0
        while (i<4):
            self.grid_rowconfigure(i, weight=1)
            i+=1
################################################################################

class waiting_window(tkinter.Frame):
    """
    It is the window where the proper analysis is taking place. The user hits
    the button and the analysis follows.

    Methods:
        analysis -> It is the method where all the computations are located. It
        starts when the user hit the button analysis.

    Attributes:

    """
    def __init__(self, parent, controller):
        """Constructor"""
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        ### title
        title = tkinter.Label(self, text="Waiting Room.",
                      font=TITLE_FONT)
        title.pack(side=tkinter.TOP)
        subtitle = tkinter.Label(self, text="""
        Please, hit the button Analysis to start the computations.
        Have a coffe meanwhile your data is being processed.
        """, font=text_fornt)
        subtitle.pack(side=tkinter.TOP)
        ### analysis
        analysis_button = tkinter.Button(self, text="Analysis",\
                command=self.analysis)
        analysis_button.pack()


        ### add a closing button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.pack()

    def analysis(self):
        """
        Basically it contains all the computations needed to perform the EDA.
        It contain the same steps used in the main CLI.
        It also generate the plots at the end. When the plots are generated
        move the software to the plot_window.
        """
        pdb_id = self.controller.app_data["pdbid"]
        pathname = self.controller.app_data["pathname"]
        pdbfile = self.controller.app_data["pdbfilename"]
        atom = self.controller.app_data["atom"]
        mode = self.controller.app_data["mode"]
        sys.stderr.write("the selcted mode is: {} ".format(mode))
        pdbalignedfile = str(pdb_id) + 'align.pdb'
        pdb_superimp = str(pathname) + str(pdb_id) + 'superimp.pdb'
        if not os.path.exists(str(pathname)+str(pdbfile)):
            pdbobj = pdb.PDBList()
            pdbobj.retrieve_pdb_file(pdb_id, pdir=str(pathname))
            sys.stderr.write("The structure {} have been \
retrieved.\n".format(pdb_id))
        atom_list = []
        if atom == 'CA':
            atom_list = ['CA']
        elif atom == 'Back':
            atom_list = ['N', 'CA', 'C', 'O']
        else:
            atom_list = ['N', 'CA', 'C', 'O']
        if mode == 'MD':
            pdbref = pdb.PDBList()
            ref_file = pdbref.retrieve_pdb_file(pdb_id, pdir=pathname)
            parser = pdb.PDBParser(QUIET=True)
            reference = parser.get_structure(pdb_id+'ref', ref_file)
            try:
                ED = eda.EDAnalysis(pdb_id, mode, atom_list, pathname+pdbfile,
                                    reference=reference)
            except (eda.WrongModeException, KeyError, ValueError):
                pass
        else:
            ED = eda.EDAnalysis(pdb_id, mode, atom_list, pathname+pdbfile)

        ED.superimpose_models()
        if mode == 'NMR':
            sys.stderr.write("Writting the superimposed file.\n")
            head = mdl.store_header_text(pathname+pdbfile)
            io = pdb.PDBIO()
            io.set_structure(ED.structure)
            io.save(pdb_superimp)
            mdl.merge_the_header(pdb_superimp, head, pathname+pdbalignedfile)
            os.remove(pdb_superimp)

        sys.stderr.write("Calculating means and coordinates\n")
        ED.createcordsarray()
        sys.stderr.write("Calculating covariance matrix\n")
        sys.stderr.write("Calculating eigenvalues and eigenvectors\n")
        ED.cal_cov()
        sys.stderr.write("Plotting eigenvalues\n")
        #pathplots = self.controller.app_data["pathplots"]
        n_plot = 30
        if ED.n < n_plot:
            n_plot = ED.n
        pathplots = pathname + 'plots/'
        plot = ED.plot_eig_wosv(n_plot)
        self.controller.app_data["plot"] = plot

        ### generating the new trajectory
        sys.stderr.write("Generating eigenvector trajectories\n")
        moved = ED.move_structure(2, 1, pathname, step=0.2)
        new_moved = moved[:-5]+'.pdb'
        mdl.merge_the_header(moved, head, new_moved)
        os.remove(moved)


        RMSD_plot = ED.RMSD_res_plot(4, pathplots, origin='interface')
        self.controller.app_data["RMSD_plot"] = RMSD_plot
        self.controller.show_frame("plot_window")


################################################################################

class plot_window(tkinter.Frame):
    """
    This window integrated matplotlib package in a Tkinter interface. It uses
    the previously generated plots in the analysis method. The user can hit one
    or the other button to see the eigen_plot or the RMSD_plot.

    Methods:
        eigen_plot -> It displays the plot eigenvalue vs eigenindex. It uses the
        matplotlib integration in tkinter.
        RMSD_plot -> It shows the RMSD vs Residue plot. It uses the matplotlib
        integration in tkinter.

    Attributes:

    """
    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        self.canvas = None
        eigen_button = tkinter.Button(self, text="Eigen Vectors Plot", command=self.eigen_plot)
        eigen_button.pack()

        RMSD_button = tkinter.Button(self, text="RMSD Plot", command=self.RMSD_plot)
        RMSD_button.pack()

            ### close button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.pack(side=tkinter.BOTTOM)


    def eigen_plot(self):
        """It generates a canvas to inttegrate the already generated plot"""
        if self.canvas != None:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        plot = self.controller.app_data["plot"]
        self.canvas = FigureCanvasTkAgg(plot, master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        #
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    def RMSD_plot(self):
        """It generates a canvas to inttegrate the already generated plot"""
        if self.canvas != None:
            self.canvas.get_tk_widget().destroy()
            self.toolbar.destroy()
        plot = self.controller.app_data["RMSD_plot"]
        self.canvas = FigureCanvasTkAgg(plot, master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        #
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

################################################################################







if __name__ == "__main__":

    app = EDA_app()
    app.mainloop()
