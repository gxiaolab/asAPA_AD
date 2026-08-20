#!/usr/bin/python

import argparse
import sys
import time
from time import strftime
from collections import defaultdict as dd

print(sys.version)
sys.path.append('./')
#from Sequence import *

import warnings
warnings.filterwarnings('ignore')

import os
import sys

import numpy as np
#np.random.seed(435) # for reproducibility
import pandas as pd
import pybedtools

import matplotlib.pyplot as plt
import seaborn as sns

import keras.backend as K
from keras.models import Model, load_model

import sys
##sys.path.append("scripts/")
sys.path.append("path/to/DeepRiPe/scripts")
from IntegratedGradients import *
from util_funcs import *
from plotseqlogo import seqlogo_fig_cross

import warnings
warnings.filterwarnings('ignore')

def precision(y_true, y_pred):
                true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
                #TPs=K.sum(K.round(K.clip(y_true * y_pred , 0, 1)))
                predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
                precision = true_positives / (predicted_positives + K.epsilon())
                return precision

def recall(y_true, y_pred):
                true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
                #TPs=K.sum(K.round(K.clip(y_ture * y_pred , 0, 1)))
                possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
                recall = true_positives / (possible_positives + K.epsilon())
                return recall

def encode_variant_bedline(bedline,genomefasta,flank_size=100):
        mut_a = bedline[4].split("/")[1]
        strand = bedline[5]
        if len(mut_a)==1:
                wild = pybedtools.BedTool(bedline[0] + "\t" + str(int(bedline[1])-flank_size) + "\t"  + str(int(bedline[2])+flank_size) + "\t" +
                                                                  bedline[3] + "\t" + str(mut_a) + "\t" + bedline[5], from_string=True )
                if strand == "-" :
                        mut_pos= flank_size
                else:
                        mut_pos= flank_size-1

                #wild = pybedtools.BedTool(bedline[0] + "\t" + bedline[1] + "\t" + bedline[2] + "\t" + bedline[3] + "\t"+ bedline[4] + "\t" + bedline[5], from_string=True)
                wild = wild.sequence(fi=genomefasta, tab=True, s=True)
                fastalist = open(wild.seqfn).read().split("\n")
                del fastalist[-1]
                seqs=[fasta.split("\t")[1] for fasta in fastalist]
                mut=seqs[0]
                mut = list(mut)
                mut[mut_pos] = mut_a
                mut = "".join(mut)
                seqs.append(mut)
                encoded_seqs =np.array([seq_to_1hot(seq) for seq in seqs])
                encoded_seqs = np.transpose(encoded_seqs,axes=(0,2,1))
                return(encoded_seqs)

##Function to score variants using models with both sequence and region
def score_variant_withregion(model, RBPnames, variant_bed, genomefasta, tr=0.1):
        reg_coded = np.full((250, 4), 0.25)
        region=np.array([reg_coded,reg_coded])
        ref_score_list = []
        alt_score_list = []
        rbp_list = []
        bed_list = []
        for bedline in variant_bed:
                #print("Bedline:", bedline)
                chrm, pos1, pos2, snp_label, ref_alt, strd = str(bedline).strip().split()
                snp_pos = '{}_{}_{}'.format(chrm, pos1, strd)
                for rbp in RBPnames:
                        RBP_index = np.where(RBPnames == rbp)[0][0]
                        encoded_seqs = encode_variant_bedline(bedline,genomefasta)
                        #print(rbp, RBP_index, encoded_seqs)
                        if encoded_seqs is not None:
                                pred = model.predict([encoded_seqs,region])[:,RBP_index]
                                #score = (pred[1]-pred[0])
                                #print(pred[1], pred[0])
                                #if min(abs(pred[0]), abs(pred[1])) >= tr: #Requiring a minimum between the 2 scores; Makes sure there is definitely binding
                                if abs(pred[0]) >= tr or abs(pred[1]) >= tr: # Only requiring one allele to have the minimum score; Catches cases where allele disrupts binding
                                        ref_score_list.append(pred[0])
                                        alt_score_list.append(pred[1])
                                        bed_list.append(bedline)
                                        rbp_list.append(rbp)
                                        #print(rbp, model, bedline, pred[0], pred[1])
        return(ref_score_list, alt_score_list, bed_list, rbp_list)

####Function to plot attribution map for variant using models with both sequence and region
def plot_variant_map_withregion(bedlines, igres, RBPnames, genomefasta, figsave_path=""):
        for bedline in bedlines[0:bedlines.count()]:
                chrm, pos1, pos2, snp_label, ref_alt, strd = str(bedline).strip().split()
                snp_pos = '{}_{}_{}'.format(chrm, pos1, strd)
                for rbp in RBPnames:
                        RBP_index = np.where(RBPnames == rbp)[0][0]
                        encoded_seqs = encode_variant_bedline(bedline,genomefasta)
                        reg_coded = np.full((250, 4), 0.25)
                        ex_seq = np.array([igres.explain([encoded_seqs[i],reg_coded],outc=RBP_index,reference=False)[0] for i in [0,1]])
                        plt.close("all")
                        seqlogo_fig_cross(np.transpose(ex_seq[:,75:125,:4],axes=(1,2,0)), vocab="RNA", figsize=(4,1.5), ncol=1, crosssite=True, cross_positions=[(24.5,25.5),(24.5,25.5),(24.5,25.5)])
                        plt.savefig(os.path.join(figsave_path, rbp + "_" + bedline[3]+".pdf"),dpi=300)
                        #plt.savefig(os.path.join(figsave_path,bedline[3]+"_"+str(RBP_index)+"_region.png"),dpi=300)
                        #plt.show()

def main(argv):
        print('job starts', strftime('%a, %d %b %Y %I:%M:%S'))
        start_time = time.time()

        genomefasta = opts.g
        variant_bed = open(opts.v).read()
        variant_bed = pybedtools.BedTool(variant_bed, from_string=True)

        model_dict = dd(list)
        with open(opts.m) as f:
                for line in f:
                        model_name, path, rbps, input_size = line.strip().split('\t')
                        #if model_name.upper().startswith("PARCLIP"): continue
                        model_dict[model_name] = [path, rbps, int(input_size)]

        with open(opts.o, 'w') as out:
                out.write("RBP\tModel\tVariant\tRef_score\tAlt_Score\n")
                for model_name in model_dict.keys():
                        model = load_model(model_dict[model_name][0], custom_objects={'precision': precision,'recall': recall })
                        #igres = integrated_gradients(model)

                        ## DEBUGGING ##
                        num_model_outputs = model.output_shape[1]
                        print(f"Model '{model_name}' has {num_model_outputs} output neurons (i.e., RBP predictions).")
                        ## DEBUGGING ##

                        RBPnames = np.array(model_dict[model_name][1].split('.'))

                        ## DEBUGGING ##
                        print(f"RBP names specified in tmp_models.txt: {RBPnames}")
                        print(f"Number of RBP names specified: {len(RBPnames)}")
                        ## DEBUGGING ##

                        ref_score_list, alt_score_list, bed_list, rbp_list = score_variant_withregion(model, RBPnames, variant_bed, genomefasta, tr=float(opts.t))
                        #print(len(ref_score_list), len(alt_score_list), len(bed_list), len(rbp_list))
                        print(ref_score_list, alt_score_list, len(bed_list), len(rbp_list))
                        for ref_score, alt_score, var, rbp in zip(ref_score_list, alt_score_list, bed_list, rbp_list):
                                out.write('{}\t{}\t{}\t{}\t{}\n'.format(rbp, model_name, str(var).strip('\n').replace('\t', '_'), str(ref_score), str(alt_score)))
                        #print(model_name, bed_list)
                        #plot_variant_map_withregion(variant_bed, igres, RBPnames, genomefasta, figsave_path=opts.p)


        print("--- %s seconds ---" % (time.time() - start_time))
        print('DONE!', strftime('%a, %d %b %Y %I:%M:%S'))

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Script descriptions here')
        parser.add_argument('-v', metavar='variant_bed', required=True, help='Variant BED')
        parser.add_argument('-m', metavar='models', required=True, help='Models')
        parser.add_argument('-g', metavar='genomefasta', required=True, help='Genome Fasta')
        parser.add_argument('-p', metavar='figsave_path', required=True, help='FigSave-Path')
        parser.add_argument('-t', metavar='threshold', required=True, help='Score threshold')
        parser.add_argument('-o', metavar='outf', required=True, help='Score distribution')

        opts = parser.parse_args()
        print('Variant-Bed: %s' % opts.v)
        print('Models: %s' % opts.m)
        print('Genome-Fasta: %s' % opts.g)
        print('FigSave-Path: %s' % opts.p)
        print('Score-Thresh: %s' % opts.t)
        print('Outf: %s' % opts.o)

        main(sys.argv[1:])
