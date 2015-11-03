
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
                #if str(motif) in line:
                    name = line.split(" ")[1]
                    flag += 1
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                #print line
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    #print words[0],
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    check += 1
        return areapwm, name


def get_chip(ids):
    chip = ChipSeq.objects.get(tf_id=ids).chipdata_set.all()
    tf_list = []
    for i in chip:
        tf_list.append(i.at_100)
    return tf_list


def get_tf(tf):
    tf_class = Matrix.objects.filter(motif_name=tf)[0].tf_id_id
    motifs = Matrix.objects.filter(tf_id=tf_class)
    #tf_class_id = motifs[0].tf_id_id
    return tf_class, motifs


