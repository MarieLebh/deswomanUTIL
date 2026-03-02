import customtkinter as ctk
from tkinter import filedialog
import subprocess
import sys
from tkinter import *
from tkinter.filedialog import *
from tkinter import filedialog
import os #so that its runnable from wherever (hopefully)

"""
This script creates a GUI for the deswomanUTIL de novo filter
"""

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("green")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("deswomanUTIL: Filter for de novo emerged neORFs GUI")
        self.geometry("1200x1000")
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        #self.img = PhotoImage(file=os.path.join(BASE_DIR,  "FilteringWorkflow.png"))  # keep reference
        #self.canvas = Canvas(self)
        #self.canvas.configure(
        #    width=self.img.width(),
        #    height=self.img.height(),
        #    bg="#FFFFFF",
        #    highlightbackground="#676869",
        #)
        #self.canvas.create_image(self.img.width()/2 , self.img.height()/2, image=self.img)
        #self.canvas.place(relx=0.2, rely=0.85, anchor="center")


        #Lets define the scripts inputs
        self.inputs = {}

        #Define the Argparse input
        self.input_defs = [
            #Basic parameters
            ("DESwoMAN Output Folder*", "--deswoman", "folder"),
            ("Path of Results folder (Filtered neORFs)*", "--out", "text"),
            ("Species tree (.nwk)*", "--tree", "file"),
            ("Sample IDs (.txt)*", "--species_file", "file"),
            ("Accepted mutations (comma separated list)*", "--accepted_mutations", "text"),
            ("Frameshift score*", "--frameshift_score", "text"),
            ("Orthofinder file (.txt)", "--ortho", "file"),
            ("Filtering type (which sequence type to use for BLAST)", "--filter_type", "choice", ["neORF", "transcript"]),
            #All TE stuff
            ("\n\nTE BLAST?\n\n", "--te_check", "bool"),
            ("TE database (.fasta)", "--te_db", "file"),
            ("Evalue (TEBlast)", "--te_eval", "choice", ["0.05", "0.01", "0.001", "0.0001"]),
            ("Coverage Blast TE (0 - 100 %)", "--te_cov", "text"),
            ("Percent Identity Blast TE (0 - 100 %)", "--te_idt", "text"),
            #All TE stuff
            ("\n\nTRANSCRIPT BLAST\n\n", "--rna_check", "bool"),
            ("Transcript (cDNA/ncRNA) database (.fasta)", "--tr_db", "file"),
            ("Evalue (Transcript Blast)", "--tr_eval", "choice", ["0.05", "0.01", "0.001", "0.0001"]),
            ("Coverage Blast Transcripts (0 - 100 %)", "--tr_cov", "text"),
            ("Percent Identity Blast Transcripts (0 - 100 %)", "--tr_idt", "text"),
            ("Transcript Blast Direction", "--tr_strand", "choice", ["plus", "both", "minus"])
        ]

        #Define the scripts (i.e. which conversion do you want!)
        self.script_map = {
            "Filter the DESwoMAN output": os.path.join(BASE_DIR,"filterDESwoMAN.py")
        }

        #Now thefine which arguments the respective script accepts
        self.script_args = {
            "Filter the DESwoMAN output": {
                "--deswoman", "--te_db", "--ortho", "--filter_type", "--te_cov", "--te_eval", "--te_idt", "--tr_cov", "--tr_idt", "--tr_eval","--tr_db", "--tr_strand", "--rna_check", "--te_check", "--out", "--tree", "--species_file" , "--accepted_mutations", "--frameshift_score"
            }
        }

        self.build_inputs()
        self.build_script_selector()
        self.build_run_button()
        self.update_visible_inputs()


    #Define the input
    def build_inputs(self):
        for item in self.input_defs:
            self.add_input_row(*item)

    def add_input_row(self, label, arg, kind, choices=None):
        frame = ctk.CTkFrame(self)
        frame.pack(fill="x", padx=20, pady=5)

        ctk.CTkLabel(frame, text=label, width=200, anchor="w").pack(side="left")

        if kind == "bool":
            var = ctk.BooleanVar(value=False)
        else:
            var = ctk.StringVar()

        self.inputs[arg] = {"var": var, "frame": frame}

        if kind == "choice":
            var.set(choices[0])
            ctk.CTkOptionMenu(
                frame,
                values=choices,
                variable=var
            ).pack(side="left", padx=10)

        elif kind == "text":
            ctk.CTkEntry(
                frame,
                textvariable=var,
                placeholder_text="Enter value..."
            ).pack(side="left", fill="x", expand=True, padx=10)

        elif kind == "bool":
            ctk.CTkCheckBox(
                frame,
                text=label,
                variable=var
            ).pack(side="left", padx=10)

        else:  # file
            ctk.CTkEntry(
                frame,
                textvariable=var
            ).pack(side="left", fill="x", expand=True, padx=10)

            ctk.CTkButton(
                frame,
                text="Select Path",
                width=80,
                command=lambda v=var, k=kind: self.browse(v, k)
            ).pack(side="right")

    #Create the script selector
    def build_script_selector(self):
        frame = ctk.CTkFrame(self)
        frame.pack(pady=20)

        self.script_var = ctk.StringVar(
            value=list(self.script_map.keys())[0]
        )

        self.script_var.trace_add("write", lambda *_: self.update_visible_inputs())

        ctk.CTkLabel(frame, text="Operation to run:").pack(side="left", padx=10)
        ctk.CTkOptionMenu(
            frame,
            values=list(self.script_map.keys()),
            variable=self.script_var
        ).pack(side="left")

    #Add run button
    def build_run_button(self):
        self.run_button = ctk.CTkButton(
            self,
            text="Filter your data!",
            height=60,
            command=self.run_script
        )
        self.run_button.pack(pady=20)

    #Define browsing

    def browse(self, var, kind):
        if kind == "file":
            path = filedialog.askopenfilename()
        elif kind == "folder":
            path = filedialog.askdirectory()
        else:
            return

        if path:
            var.set(path)

    def update_visible_inputs(self):
        script = self.script_var.get()
        allowed_args = self.script_args.get(script, set())

        for arg, data in self.inputs.items():
            frame = data["frame"]
            var = data["var"]

            if arg in allowed_args:
                frame.pack(fill="x", padx=20, pady=5)
            else:
                frame.pack_forget()
                if isinstance(var, ctk.BooleanVar):
                    var.set(False)
                else:
                    var.set("")

    def run_script(self):
        script_name = self.script_var.get()
        script_path = self.script_map[script_name]

        cmd = [sys.executable, script_path]

        allowed_args = self.script_args.get(script_name, set())

        for arg in allowed_args:
            data = self.inputs[arg]
            value = data["var"].get()

            if isinstance(data["var"], ctk.BooleanVar):
                if value:
                    cmd.append(arg)
            else:
                if value:
                    cmd.extend([arg, value])


        print("Running:", cmd)
        subprocess.Popen(cmd)
        self.destroy()  


if __name__ == "__main__":
    app = App()
    app.mainloop()
