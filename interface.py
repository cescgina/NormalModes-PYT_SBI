import tkinter
import module_david as md
import os.path
import sys
from tkinter import filedialog
import Bio.PDB as pdb

class initial_root:
    def __init__(self, master):
        self.master = master
        ### main window
        master.title("EDA: by JF Gilabert & D Mas")
        master.geometry("600x400")
        ### title
        self.title = tkinter.Label(master, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25))
        self.title.pack(side=tkinter.TOP)
        ### photo
        self.photo = tkinter.PhotoImage(file="image.gif")
        self.w = tkinter.Label(master, image=self.photo)
        self.w.photo = self.photo
        self.w.pack(side=tkinter.TOP)
        ### subtitle
        self.subtitle = tkinter.Label(master, text="Please, enter a file or a pdb code.",
                                 font=("Helvetica", 20))
        self.subtitle.pack(side=tkinter.TOP)
        ### pdb code entry
        self.labelframe_code = tkinter.Frame(master, bd=2)
        self.labelframe_code.pack(side=tkinter.TOP)
        self.entry_var = tkinter.StringVar()
            ### entry
        self.label = tkinter.Label(self.labelframe_code, text="Enter a pdb code:",
                              font=("Helvetica", 15), width=40)
        self.label.pack(side=tkinter.TOP)
        self.entry = tkinter.Entry(self.labelframe_code, bd=2, width=8, textvariable=self.entry_var)
        self.entry.pack(side=tkinter.TOP)
            ### check button
        self.b = tkinter.Button(self.labelframe_code, text="Check",
                           command=self.check_uniprot_accession_code)
        self.b.pack(side=tkinter.TOP)
            ### get pdb button
        self.get = tkinter.Button(self.labelframe_code, text="Get PDB", command=self.get_pdb)
        self.get.pack(side=tkinter.TOP)

            ### add a file button
        self.add = tkinter.Button(self.labelframe_code, text="Add a pdb file", command=self.select_file)
        self.add.pack(side=tkinter.TOP)
            ### close button
        self.close_button = tkinter.Button(master, text="Close", command=master.quit)
        self.close_button.pack()

    def check_uniprot_accession_code(self):
        if md.pdb_code_check(self.entry_var.get()):
            self.entry["foreground"] = "green"
        else:
            sys.stderr.write("%s is not a uniprotaccession code format\n" \
                            %self.entry_var.get())
            self.entry["foreground"] = "red"
    def get_pdb(self):
        #structure = md.from_pdb_code_to_structure(str(self.entry_var.get()))
        #sys.stderr.write("The structure {} have been retrieved.\n".format(self.entry_var.get()))
        global interface_code
        interface_code = str(self.entry_var.get())
        self.master.destroy()


    def select_file(self):
        filename = filedialog.askopenfilename()
        sys.stderr.write("Reading the structure from {}.\n".format(filename))
        global path_name, pdb_file
        path_name,pdb_file = os.path.split(filename)
        self.master.destroy()

class waiting_window:
    """docstring for waiting_window"""
    def __init__(self, master):
        self.master = master
        ### main window
        master.title("EDA: by JF Gilabert & D Mas")
        master.geometry("600x400")
        ### title
        self.title = tkinter.Label(master, text="Waiting Room.",
                      font=("Helvetica", 25))
        self.title.pack(side=tkinter.TOP)
        self.subtitle = tkinter.Label(master, text="Please, have a coffe meanwhile \
your data is being processed.", font=("Helvetica", 20))
        self.subtitle.pack(side=tkinter.TOP)

        self.close_button = tkinter.Button(master, text="Close", command=master.quit)
        self.close_button.pack()




if __name__ == "__main__":
    root = tkinter.Tk()
    my_gui = initial_root(root)
    root.mainloop()
    if interface_code:
        structure = md.from_pdb_code_to_structure(interface_code)
        sys.stderr.write("The structure {} have been retrieved.\n".format(interface_code))
        sys.stderr.write("The structure {} have been loaded \n".format(structure.header['name']))
        root=tkinter.Tk()
        waiting_gui = waiting_window(root)
        root.mainloop()


    elif pdb_file:
        parser = pdb.PDBParser(QUIET=True)
        pdb_structure = parser.get_structure(pdb_file, pdb_file)
        sys.stderr.write("The file {} have been written\n".format((pdb_file)))
    else:
        sys.stderr.write("Neither file or code provided\n")
