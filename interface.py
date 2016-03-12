import tkinter
import module_david as md
import sys
from tkinter import filedialog

# PLACING THE ROOT WINDOW
root = tkinter.Tk()

# SETTING THE width and the wide
root.geometry("600x400")
# Setting the title
root.wm_title("EDA: by JF Gilabert & D Mas")

# Generate a title for the WINDOW

title = tkinter.Label(root, text="Welcome to the ED Analyser.",
                      font=("Helvetica", 25))
title.pack(side=tkinter.TOP)

photo = tkinter.PhotoImage(file="image.gif")
w = tkinter.Label(root, image=photo)
w.photo = photo
w.pack(side=tkinter.TOP)

subtitle = tkinter.Label(root, text="Please, enter a file or a pdb code.",
                         font=("Helvetica", 20))
subtitle.pack(side=tkinter.TOP)


# Generate a pdb code entry button.

labelframe_code = tkinter.Frame(root, bd=2)
labelframe_code.pack(side=tkinter.TOP)

entry_var = tkinter.StringVar()


def check_uniprot_accession_code():
    global entry_var
    global entry
    if md.pdb_code_check(entry_var.get()):
        entry["foreground"] = "green"
    else:
        sys.stderr.write("%s is not a uniprotaccession code format\n" %
                         entry_var.get())
        entry["foreground"] = "red"
# def print_env():
#     global entry_var
#     global entry
#     sys.stderr.write("%s\n" %entry_var.get())


def get_pdb():
    global entry_var
    global entry
    md.from_pdb_code_to_structure(str(entry_var.get()))
    print("Finnished!\n")


def quit_window():
    sys.exit()

# CREATE THE LABEL
label = tkinter.Label(labelframe_code, text="Enter a pdb code:",
                      font=("Helvetica", 15), width=40)
label.pack(side=tkinter.TOP)

# CREATE THE ENTRY
entry = tkinter.Entry(labelframe_code, bd=2, width=8, textvariable=entry_var)
entry.pack(side=tkinter.TOP)

b = tkinter.Button(labelframe_code, text="Check",
                   command=check_uniprot_accession_code)
b.pack(side=tkinter.TOP)

s = tkinter.Button(labelframe_code, text="Submit", command=get_pdb)
s.pack(side=tkinter.TOP)

# GETTING THE file


def select_file():
    filename = filedialog.askopenfilename()

s = tkinter.Button(labelframe_code, text="Add a pdb file", command=select_file)
s.pack(side=tkinter.TOP)

# GENERATING THE MAIN LOOP
root.mainloop()
