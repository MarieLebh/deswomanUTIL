import customtkinter as ctk
from tkinter import filedialog
import subprocess
import sys
from tkinter import *
from tkinter.filedialog import *
from tkinter import filedialog
import os #so that its runnable from wherever (hopefully)

"""
This script creates a GUI for the deswomanUTIL file conversion scripts!
"""

ctk.set_appearance_mode("System")
ctk.set_default_color_theme("green")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("deswomanUTIL: File Converter GUI")
        self.geometry("800x500")
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        #self.img = PhotoImage(file=os.path.join(BASE_DIR,  "Logo.png"))  # keep reference
        #self.canvas = Canvas(self)
        #self.canvas.configure(
        #    width=self.img.width(),
        #    height=self.img.height(),
        #    bg="#FFFFFF",
        #    highlightbackground="#676869",
        #)
        #self.canvas.create_image(self.img.width()/2 , self.img.height()/2, image=self.img)
        #self.canvas.place(relx=0.5, rely=0.6, anchor="center")


        #Lets define the scripts inputs
        self.inputs = {}

        #Define the Argparse input
        self.input_defs = [
            ("DESwoMAN Information File (.txt)", "--deswoman", "file"),
            ("Transcriptome Assembly File (.gtf)", "--gtf", "file"),
            ("Choice (feature to keep)", "--choice", "choice", ["transcript", "neORF", "both"]),
            ("Output Name (without extension)", "--outname", "text"),
            ("Add stringtie locus (based on gene id)", "--add_stringtie_locus", "bool"), 
            ("Collaps de novo ORFs", "--collapse_orf", "bool"),
            ("Collapse de novo transcripts", "--collapse_denovo", "bool") 
        ]

        #Define the scripts (i.e. which conversion do you want!)
        self.script_map = {
            "TransformInfoFiletoBed": os.path.join(BASE_DIR,"getBEDfromOut.py"),
            "TransformInfoFiletoGFF": os.path.join(BASE_DIR,  "getGFFfromOut.py"),
        }

        #Now thefine which arguments the respective script accepts
        self.script_args = {
            "TransformInfoFiletoBed": {
                "--deswoman", "--gtf", "--choice", "--outname"
            },
            "TransformInfoFiletoGFF": {
                "--deswoman", "--gtf", "--outname", "--add_stringtie_locus", "--collapse_orf", "--collapse_denovo"
            },
        }

        self.build_inputs()
        self.build_script_selector()
        self.build_run_button()
        self.update_visible_inputs()

    def build_inputs(self):
        """
        Defines the input
        """
        for item in self.input_defs:
            self.add_input_row(*item)

    def add_input_row(self, label, arg, kind, choices=None):
        """
        Defines the input kinds
        """
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
                text="Select File",
                width=80,
                command=lambda v=var, k=kind: self.browse(v, k)
            ).pack(side="right")

    def build_script_selector(self):
        """
        Create the script selector (i.e. bed vs gff transformation as done here)
        """
        frame = ctk.CTkFrame(self)
        frame.pack(pady=20)

        self.script_var = ctk.StringVar(
            value=list(self.script_map.keys())[0]
        )

        self.script_var.trace_add("write", lambda *_: self.update_visible_inputs())

        ctk.CTkLabel(frame, text="Select File Conversion:").pack(side="left", padx=10)
        ctk.CTkOptionMenu(
            frame,
            values=list(self.script_map.keys()),
            variable=self.script_var
        ).pack(side="left")

    def build_run_button(self):
        """
        Add run button and defines what it sais
        """
        self.run_button = ctk.CTkButton(
            self,
            text="Start File Conversion",
            height=60,
            command=self.run_script
        )
        self.run_button.pack(pady=20)

    def browse(self, var, kind):
        """
        Define browsing
        """
        if kind == "file":
            path = filedialog.askopenfilename()
        elif kind == "folder":
            path = filedialog.askdirectory()
        else:
            return

        if path:
            var.set(path)

    def update_visible_inputs(self):
        """
        Modify which input is visible based on script to run
        """
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
        """
        Run the script and exit when done
        """
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
