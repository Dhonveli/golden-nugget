#!/usr/bin/env python


import argparse
import os 
import sys


MATCH = +3
MISMATCH = -3
GAP = -2
def params():
    obj_parse_args = argparse.ArgumentParser(description="""Given two sequences or a file containing sequences return the alignments. 
        The sequences in the input file have to be tab separated""")
    obj_parse_args.add_argument("--seq1",type=str,help="The first sequence to align")
    obj_parse_args.add_argument("--seq2",type=str,help="The second sequence to align")
    obj_parse_args.add_argument("-i","--input",type=str,help="The file in input")
    obj_parse_args.add_argument("-m","--match",default=MATCH,type=float,required=False,
                help="score of the match (default = {})".format(MATCH))
    obj_parse_args.add_argument("-s","--mismatch",default=MISMATCH,type=float,required=False,
                help="score of the mismatch (default = {})".format(MISMATCH))
    obj_parse_args.add_argument("-g","--gap",default=GAP,type=float,required=False,
                help="score of the gap (default = {})".format(GAP))
    obj_parse_args.add_argument("--minscore",type=float,help="The minimum score")
    obj_parse_args.add_argument("--minlength",type=int,help="The minimum length")
    obj_parse_args.add_argument("--numresult",type=int,help="The number of alignments to return")
    args = obj_parse_args.parse_args()
    return args

def check_params(parameters):
    parameters.seqs = []
    input_f = False
    if parameters.input is not None:
        if os.path.exists(parameters.input):
            if os.path.isfile(parameters.input):
                size = len(parameters.seqs)
                parameters.seqs += [tuple(line.strip().split("\t")) for line 
                    in open(parameters.input) if len(line.strip().split("\t"))==2]
                input_f = True
                if len(parameters.seqs) <= size:
                    print("""the file could be empty or wrongly formatted: \nseq1\tseq2
                        \nseq3\tseq4\n...\n(tab separated)""")
            else:
                print("input is not a file")
        else:
            print("input file doesn't exist")
    if ((parameters.seq1 is not None) and (parameters.seq2 is not None)):
            parameters.seqs.append((parameters.seq1,parameters.seq2))
    if ((parameters.seq1 is None) or  (parameters.seq2 is None)) and input_f is False:
        print("You have to give in input at least two sequences or a file of sequences !")            
        if not parameters.seqs:
            print("You did not give me anything...")
            sys.exit(1)
    if parameters.match <= 0:
        print("provided value for matches is not correct, default will be used instead")
        parameters.match = MATCH
    if parameters.mismatch > 0:
        print("provided value for mismatches is not correct, default will be used instead")
        parameters.mismatch= MISMATCH
    if parameters.gap > 0:
        print("provided value for gaps is not correct, default will be used instead")
        parameters.gap = GAP
    if parameters.minscore is not None:
        if parameters.minscore < 0:
            print("provided value for minimium score is negative, no minimum score will be used")
            parameters.minscore = None
    if parameters.minlength is not None:
        if parameters.minlength < 0:
            print("provided value for minimium length is negative, no minimum length will be used")
            parameters.minlength = None
    return parameters
    

            
def smith_waterman(seqs,m=MATCH,ms=MISMATCH,g=GAP):
    seq1 = seqs[0]
    seq2 = seqs[1]
    n = len(seqs[0])
    t = len(seqs[1])
    DP = [[0 for x in range(n+1)] for y in range(t+1)] # matrix of value
    DM = [[[0,0,0] for x in range(n+1)] for y in range(t+1)] # matrix of parameters # D U L
    dic_sc = {}
    for i in range(t+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                DP[i][j] == 0
            else: 
                if seq2[i-1] == seq1[j-1]:
                    DP[i][j] = max(DP[i-1][j-1] + m, DP[i-1][j] + g, DP[i][j-1] + g,0)   
                    if DP[i][j] == DP[i-1][j-1] + m: # match
                        DM[i][j][0] = 1
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]
                    if DP[i][j] == DP[i-1][j] + g: # gap
                        DM[i][j][1] = 1
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]
                    if DP[i][j] == DP[i][j-1] + g: # gap
                        DM[i][j][2] = 1
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)] 
                    elif DP[i][j] == 0:         
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]
                else:
                    DP[i][j] = max(DP[i-1][j-1] + ms, DP[i-1][j] + g, DP[i][j-1] + g,0)
                    if DP[i][j] == DP[i-1][j-1] + ms: # mismatch
                        DM[i][j][0] = 1 
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]                    
                    if DP[i][j] == DP[i-1][j] + g: # gap
                        DM[i][j][1] = 1   
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]              
                    if DP[i][j] == DP[i][j-1] + g: # gap
                        DM[i][j][2] = 1
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]
                    if DP[i][j] == 0:
                        if DP[i][j] in dic_sc:
                            dic_sc[DP[i][j]].append((i,j))
                        else:
                            dic_sc[DP[i][j]] = [(i,j)]
    return(DP,DM,dic_sc)

def print_align(DP,DM,seqs,coord):
    seq1 = seqs[0]
    seq2 = seqs[1]
    seq_list = []
    score = 1
    tmp_j = coord[0]
    tmp_i = coord[1]
    list3 = []
    list1 = []
    list2 = []
    while score != 0:
        if DM[tmp_j][tmp_i][1] == 1:
            tmp_j = tmp_j-1
            score = DP[tmp_j][tmp_i]
            list1.append("_")
            list2.append(seq2[tmp_j])
            list3.append(" ")
        elif DM[tmp_j][tmp_i][2] == 1:
            tmp_i = tmp_i-1
            score = DP[tmp_j][tmp_i]    
            list2.append("_")   
            list1.append(seq1[tmp_i])
            list3.append(" ")
        elif DM[tmp_j][tmp_i][0] == 1:
            tmp_i = tmp_i-1
            tmp_j = tmp_j-1
            score = DP[tmp_j][tmp_i]
            list2.append(seq2[tmp_j])
            list1.append(seq1[tmp_i])
            if seq1[tmp_i] == seq2[tmp_j]:
                list3.append("|")
            else:
                list3.append(":")
    seq_list.append((list1,list2,list3))  
    return seq_list


if __name__ == "__main__":
    args = params()
    args = check_params(args)
    for seqs in args.seqs:
        DP,DM,dic_sc = smith_waterman(seqs,args.match,args.mismatch,args.gap)
        scores = sorted(list(dic_sc.keys()),reverse = True)
        tmp = len(dic_sc[scores[0]])
        if args.numresult is not None:
            count = args.numresult
            for i in scores:
                if args.minscore is not None:
                    if args.minscore < i:
                        if count > 0: 
                            for j in dic_sc[i]:
                                seq_list = print_align(DP,DM,seqs,j)
                                for seq_tu in seq_list:
                                    seq1 = "".join(seq_tu[0][::-1])
                                    seq2 = "".join(seq_tu[1][::-1])
                                    mid = "".join(seq_tu[2][::-1])
                                    if args.minlength is not None:
                                        if args.minlength < len(seq1):
                                            count -= 1
                                            print("SCORE: ",i)
                                            print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))

                                    else:
                                        count -= 1
                                        print("SCORE", i)
                                        print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))
                else:
                    if count > 0:
                        for j in dic_sc[i]:
                            seq_list = print_align(DP,DM,seqs,j)
                            for seq_tu in seq_list:
                                    seq1 = "".join(seq_tu[0][::-1])
                                    seq2 = "".join(seq_tu[1][::-1])
                                    mid = "".join(seq_tu[2][::-1])
                                    if args.minlength is not None:
                                        if args.minlength < len(seq1):
                                            count -= 1
                                            print("SCORE: ",i)
                                            print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))

                                    else:
                                        count -= 1
                                        print("SCORE", i)
                                        print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))
        else:
            for i in scores:
                if args.minscore is not None:
                    if args.minscore < i:
                        for j in dic_sc[i]:
                            seq_list= print_align(DP,DM,seqs,j)
                            for seq_tu in seq_list:
                                seq1 = "".join(seq_tu[0][::-1])
                                seq2 = "".join(seq_tu[1][::-1])
                                mid = "".join(seq_tu[2][::-1])
                                if args.minlength is not None:
                                    if args.minlength < len(seq1):
                                        count -= 1
                                        print("SCORE: ",i)
                                        print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))

                                else:
                                    count -= 1
                                    print("SCORE", i)
                                    print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))
                else:
                    for j in dic_sc[i]:
                        seq_list = print_align(DP,DM,seqs,j)
                        for seq_tu in seq_list:
                            seq1 = "".join(seq_tu[0][::-1])
                            seq2 = "".join(seq_tu[1][::-1])
                            mid = "".join(seq_tu[2][::-1])
                            if args.minlength is not None:
                                if args.minlength < len(seq1):
                                    print("SCORE: ",i)
                                    print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))

                            else:
                                print("SCORE", i)
                                print("\n {} \n {} \n {} \n ".format(seq1,mid,seq2))
                break