'''
This is a module to axtract a motif from meme file given a unique motif name 
and create dictionary for sequence scoring 

Usage:
    python fetch_motif.py <meme_file> <Motfi_name (optional)>

'''
import sys


def get_motif_from_meme(meme, motif="MOTIF"):
    """
    Extract a motif from meme file given a unique motif
    name and create dictionary for sequence scoring

    Default motif name is keyword MOTIF for single motif files.
    """
    name = ""
    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    flag = 0
    check = 0
    with open(meme, "r") as f1:
        for line in f1:
            if line.startswith('MOTIF'):
                if line.split(" ")[1] == motif:
                    name = line.split(" ")[1]
                    flag += 1
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    check += 1
        return areapwm, name

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(1)
    meme_file = sys.argv[1]
    motif_name = sys.argv[2]
    get_motif_from_meme(meme, motif="MOTIF")

