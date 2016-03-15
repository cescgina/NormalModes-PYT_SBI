import tkinter
import module_david as mdl
import os.path
import os
import sys
from tkinter import filedialog
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import Bio.PDB as pdb
import edanalysis as eda

# implement the default mpl key bindings


TITLE_FONT = ("Helvetica", 18, "bold")
text_fornt = ("Courier", 14)
###############################################################################
class EDA_app(tkinter.Tk):

    def __init__(self, *args, **kwargs):
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
        self.title("EDA: by JF Gilabert & D Mas")
        self.geometry("600x600")
        container = tkinter.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in (StartPage, initial_root, waiting_window, About_EDA, plot_window):
            page_name = F.__name__
            frame = F(container, self)
            self.frames[page_name] = frame
            # put all of the pages in the same location;
            # the one on the top of the stacking order
            # will be the one that is visible.
            frame.grid(row=0, column=0, sticky="nsew")
        ### Adding the plot page.

        self.show_frame("StartPage")

    def show_frame(self, page_name):
        '''Show a frame for the given page name'''
        frame = self.frames[page_name]
        frame.tkraise()

###############################################################################
class StartPage(tkinter.Frame):

    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        title = tkinter.Label(self, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25))
        title.pack(side=tkinter.TOP)
        photo = tkinter.PhotoImage(file="image.gif")
        w = tkinter.Label(self, image=photo)
        w.photo = photo
        w.pack(side=tkinter.TOP)
        subtitle = tkinter.Label(self, text="This Software is the result for a \
project of Structual Bioinformatics \n and Introduction to Python from the Master \
of Science in Bioinformatics \n at Universitat Pompeu Fabra. We use the Essential \
Dynamics approach \n in order to analyse the Normal Modes of NMR structures and \
Molecular \n Dynamics. The input format is pdb although it  can be \n also retrieved \
automatically by the software.", font=text_fornt)
        subtitle2 = tkinter.Label(self, text="Authors: JF Gilabert & D Mas",\
         font=("Helvetica", 22, "bold"))
        subtitle.pack(side=tkinter.TOP)
        subtitle2.pack(side=tkinter.TOP)

        button1 = tkinter.Button(self, text="Start the Analysis",
                            command=lambda: controller.show_frame("initial_root"))
        button2 = tkinter.Button(self, text="About EDA",
                            command=lambda: controller.show_frame("About_EDA"))
        button1.pack()
        button2.pack()



################################################################################

class initial_root(tkinter.Frame):
    """
    This is the initial page of the application where the files are loaded
    """
    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        ### title
        title = tkinter.Label(self, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25))
        title.pack(side=tkinter.TOP)
        ### photo
        photo = tkinter.PhotoImage(file="image.gif")
        w = tkinter.Label(self, image=photo)
        w.photo = photo
        w.pack(side=tkinter.TOP)
        ### subtitle
        subtitle = tkinter.Label(self, text="Please, enter a file or a pdb code.",
                                 font=("Helvetica", 20))
        subtitle.pack(side=tkinter.TOP)
        ### pdb code entry
        labelframe_code = tkinter.Frame(self, bd=2)
        labelframe_code.pack(side=tkinter.TOP)
        self.entry_var = tkinter.StringVar()
            ### entry
        label = tkinter.Label(labelframe_code, text="Enter a pdb code:",
                              font=("Helvetica", 15), width=40)
        label.pack(side=tkinter.TOP)
        self.entry = tkinter.Entry(labelframe_code, bd=2, width=8, textvariable=self.entry_var)
        self.entry.pack(side=tkinter.TOP)
            ### check button
        b = tkinter.Button(labelframe_code, text="Check",
                           command=self.check_uniprot_accession_code)
        b.pack(side=tkinter.TOP)
            ### get pdb button
        get = tkinter.Button(labelframe_code, text="Get PDB", command=self.get_pdb)
        get.pack(side=tkinter.TOP)

            ### add a file button
        add = tkinter.Button(labelframe_code, text="Add a pdb file", command=self.select_file)
        add.pack(side=tkinter.TOP)
        ### choose atoms
        self.m = tkinter.StringVar()
        m1 = tkinter.Radiobutton( self, text="NMR", variable = self.m, value = 'NMR', command = self.select_m)
        m2 = tkinter.Radiobutton( self, text="MD", variable = self.m, value = 'MD', command = self.select_m )

        m1.pack( side=tkinter.TOP)
        m2.pack( side=tkinter.TOP )

        self.label_m = tkinter.Label(self)
        self.label_m.pack()
        ### choose mode
        self.var = tkinter.StringVar()
        r1 = tkinter.Radiobutton( self, text="Only CA", variable = self.var, value = 'CA', command = self.select)
        r2 = tkinter.Radiobutton( self, text="Backbone atoms", variable = self.var, value = 'Back', command = self.select )
        r3 = tkinter.Radiobutton( self, text="all", variable = self.var, value = 'all', command = self.select )

        r1.pack( side=tkinter.TOP)
        r2.pack( side=tkinter.TOP )
        r3.pack( side=tkinter.TOP )

        self.label = tkinter.Label(self)
        self.label.pack()

        ### start page button
        button = tkinter.Button(self, text="Go to the start page",
                           command=lambda: controller.show_frame("StartPage"))
        button.pack()
            ### close button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.pack()

    def select_m(self):
        selection = "You selected the option %s" %self.m.get()
        self.label_m.config(text=selection)
        self.controller.app_data["mode"] = self.m.get()

    def select(self):
        selection = "You selected the option %s" %self.var.get()
        self.label.config(text=selection)
        self.controller.app_data["atom"] = self.var.get()

    def check_uniprot_accession_code(self):
        if mdl.pdb_code_check(self.entry_var.get()):
            self.entry["foreground"] = "green"
        else:
            sys.stderr.write("%s is not a uniprotaccession code format\n" \
                            %self.entry_var.get())
            self.entry["foreground"] = "red"
    def get_pdb(self):
        interface_code = str(self.entry_var.get())
        if mdl.pdb_code_check(interface_code):
            #structure = mdl.from_pdb_code_to_structure(interface_code)
            #sys.stderr.write("The structure {} have been retrieved.\n".format(interface_code))
            self.controller.app_data["pdbid"] = interface_code
            self.controller.app_data["pdbfilename"] = 'pdb'+interface_code+ '.ent'
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

    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        label1 = tkinter.Label(self, text="Help and Documentation", font=TITLE_FONT)
        label2 = tkinter.Label(self, text="EDA is a python based software \
that performs a Normal Mode Analysis \n using a Essential Dynamics Simplification.\
", font=text_fornt)
        label3 = tkinter.Label(self, text="You can perform the analysis from a \
pdb file or you can introduce the \n desired pdb code and the program will\
retrieve the pdb file for you if \n you have internet connection.", font=text_fornt)
        label4 = tkinter.Label(self, text="The Program will output a set of files \
 and plots. First, it will \n generate a superimposed pdb files with all the NMR models\
 or MD frames. \n After that it will calculate the covariance matrix of the coordinates \n \
 and will output a plot of the eigenvalues derived from the \n diagonalitzation\
 of the matrix. Finally, it will output a pdb file\n with the moved coordinates.",\
        font=text_fornt)
        label1.pack(side="top", fill="x", pady=10)
        label2.pack(side="top", fill="x", pady=10)
        label3.pack(side="top", fill="x", pady=10)
        label4.pack(side="top", fill="x", pady=10)
        button = tkinter.Button(self, text="Go to the start page",
                           command=lambda: controller.show_frame("StartPage"))
        button.pack()
        button1 = tkinter.Button(self, text="Start the Analysis",
                            command=lambda: controller.show_frame("initial_root"))
        button1.pack()

################################################################################

class waiting_window(tkinter.Frame):
    """docstring for waiting_window"""
    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        ### title
        title = tkinter.Label(self, text="Waiting Room.",
                      font=TITLE_FONT)
        title.pack(side=tkinter.TOP)
        subtitle = tkinter.Label(self, text="Please, hit the button Analysis to \
start the computations. \n Have a coffe meanwhile your data is being \
processed.", font=text_fornt)
        subtitle.pack(side=tkinter.TOP)
        ### analysis
        analysis_button = tkinter.Button(self, text="Analysis", command=self.analysis)
        analysis_button.pack()
        ### add a closing button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.pack()

    def analysis(self):
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
        pathplots = self.controller.app_data["pathplots"]
        n_plot = 30
        if ED.n < n_plot:
            n_plot = ED.n
        pathplots = pathname + 'plots/'
        plot = ED.plot_eig_wosv(n_plot)
        self.controller.app_data["plot"] = plot


        RMSD_plot = ED.RMSD_res_plot(4, pathplots, origin='interface')
        self.controller.app_data["RMSD_plot"] = RMSD_plot
        self.controller.show_frame("plot_window")


################################################################################

class plot_window(tkinter.Frame):
    def __init__(self, parent, controller):
        tkinter.Frame.__init__(self, parent)
        self.controller = controller
        eigen_button = tkinter.Button(self, text="Eigen Vectors Plot", command=self.eigen_plot)
        eigen_button.pack()

        RMSD_button = tkinter.Button(self, text="RMSD Plot", command=self.RMSD_plot)
        RMSD_button.pack()

            ### close button
        close_button = tkinter.Button(self, text="Close", command=self.quit)
        close_button.pack()


    def eigen_plot(self):
        plot = self.controller.app_data["plot"]
        canvas = FigureCanvasTkAgg(plot, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        #
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
    def RMSD_plot(self):
        plot = self.controller.app_data["RMSD_plot"]
        canvas = FigureCanvasTkAgg(plot, master=self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        #
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
################################################################################







if __name__ == "__main__":

    app = EDA_app()
    app.mainloop()
