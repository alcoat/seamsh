import glob
import shutil
import os
try :
    shutil.rmtree("examples")
except :
    pass
os.makedirs("examples")

for filename in glob.glob("../examples/*.py") :
    inrst = True
    with open(filename, "r") as f:
        bname = os.path.basename(filename)[:-3]
        with open("examples/"+bname+".rst","w") as fout :
            fout.write("Download :download:`this example <../%s>`.\n"%filename)
            for i,line in enumerate(f) :
                line = line[:-1]
                lstrip = line.strip()
                if lstrip == "# %%" :
                    inrst = True
                    fout.write("\n")
                    continue
                if (not lstrip and inrst) or i == 0 or (inrst and lstrip[0] != '#'):
                    fout.write("\n.. code-block:: python\n  \n")
                    inrst = False
                    continue
                if inrst :
                    fout.write(line[2:]+"\n")
                else :
                    fout.write("  "+line+"\n")



