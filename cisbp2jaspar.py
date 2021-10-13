#!/usr/bin/env python

from Bio import motifs
import copy
import os
import pandas as pd
import re

def get_motifs_to_TFs(info_file):

    motif2tfs = []

    df = pd.read_csv("./data/CisBP-mouse/%s" % info_file,
        sep="\t", header=0, usecols=[3, 6, 8])
    df = df[df["TF_Status"] != "N"]
    df = df.groupby("Motif_ID")["TF_Name"].aggregate(set)
    for motif, tfs in df.iteritems():
        motif2tfs.append([motif, ";".join(sorted(tfs))])
    df = pd.DataFrame(motif2tfs, columns=["Motif_ID", "TF_Names"])
    df.to_csv("%s.tsv" % file2prefixes[info_file], sep="\t")

    return(set(df.Motif_ID.to_list()))

def get_JASPAR_bundle(info_file, allowed_pwms):

    pwms = []

    pwms_dir = "./data/CisBP-mouse/pwms_all_motifs/"

    for pwm_file in sorted(os.listdir(pwms_dir)):
        if pwm_file[:11] not in allowed_pwms:
            continue
        try:
            with open(os.path.join(pwms_dir, pwm_file)) as handle:
                pwm = motifs.read(handle, "pfm-four-columns")
                pwm.matrix_id = pwm_file[:11]
                pwm.name = pwm_file[:11]
                pwms.append(copy.deepcopy(pwm))
        except:
            pass

    with open("%s.jaspar" % file2prefixes[info_file], "w") as handle:
        handle.write(motifs.write(pwms, "jaspar"))

file2prefixes = {
    "TF_Information.txt": "CisBP-mouse",
    "TF_Information_all_motifs.txt": "CisBP-mouse.all_motifs",
    "TF_Information_all_motifs_plus.txt": "CisBP-mouse.all_motifs_plus"
}

for info_file in file2prefixes:
    get_JASPAR_bundle(info_file, get_motifs_to_TFs(info_file))