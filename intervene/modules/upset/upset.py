# coding: utf-8

"""
InterVene: a tool for intersection and visualization of multiple genomic region sets
Created on January 10, 2017
@author: <Aziz Khan>aziz.khan@ncmm.uio.no
"""
import sys
import os
import tempfile
import itertools
from intervene.modules.pairwise.pairwise import get_name
from pybedtools import BedTool, helpers
from intervene import helpers as hlp


def genomic_upset(options, label_names):
    '''
    Arguments:
        input_files -  List of BED files to to calculate the weights
        output - output path
        label_names - names of input files
    Takes a list of sets a list of the sizes of non-overlapping intersections between them 
    '''

    input_files = options.input
    output = options.output

    bed_cache = [BedTool(f) for f in options.input]


    kwargs = hlp.map_bedtools_options(options.bedtools_options)

    N = len(input_files)

    # Build a non-redundant union of all input intervals so that every
    # genomic locus is counted exactly once (merge-then-count approach).
    # This ensures symmetric counts: unlike anchoring on ones[0], the result
    # is independent of set ordering.
    union = bed_cache[0]
    for bed in bed_cache[1:]:
        union = union.cat(bed, postmerge=False)
    union = union.sort().merge()

    # Assign every union interval an integer bitmask (bit j = overlaps set j).
    # This requires exactly N bedtools intersect calls instead of O(2^N * N),
    # and produces far fewer temp files (reducing cleanup overhead).
    union_intervals = [(iv.chrom, iv.start, iv.end) for iv in union]
    coord_to_idx = {iv: i for i, iv in enumerate(union_intervals)}
    membership = [0] * len(union_intervals)

    for j, bed in enumerate(bed_cache):
        for iv in union.intersect(bed, u=True, **kwargs):
            idx = coord_to_idx[(iv.chrom, iv.start, iv.end)]
            membership[idx] |= (1 << j)

    # Count intervals per exclusive combination and build the weights dict.
    from collections import Counter, defaultdict
    counts = Counter(membership)
    weights = {}
    for mask, count in counts.items():
        if mask == 0:
            continue
        key = ''.join('1' if (mask >> j) & 1 else '0' for j in range(N))
        weights[key] = count

    # Save overlapping regions per combination if requested.
    if options.saveoverlaps:
        intervals_by_key = defaultdict(list)
        for i, mask in enumerate(membership):
            if mask == 0:
                continue
            key = ''.join('1' if (mask >> j) & 1 else '0' for j in range(N))
            intervals_by_key[key].append(union_intervals[i])

        for key, intervals in intervals_by_key.items():
            if len(intervals) >= options.overlapthresh:
                file_name = ''.join(
                    '_' + label_names[j] for j, b in enumerate(key) if b == '1'
                )
                file_name = key + file_name
                hlp.create_dir(f"{output}/sets")
                with open(f"{output}/sets/{file_name}.bed", 'w') as f:
                    for chrom, start, end in intervals:
                        f.write(f"{chrom}\t{start}\t{end}\n")

    helpers.cleanup()

    return weights

def list_upset(options, label_names):
    '''
    Arguments:
        input_files -  List of list files to calculate weights for upset plot
        output - output path
        label_names - names of input files
    Takes a list of sets a list of the sizes of non-overlapping intersections between them 
    '''
    input_files = options.input
    output = options.output
    S =[]
    for f in input_files:
        with open(f) as f_open:
            S.append(set(f_open.read().splitlines()))
    N = len(S)
    # Generate a truth table of intersections to calculate 
    truth_table = [x for x in itertools.product("01", repeat=N)][1:]
    weights = {}
    for t in truth_table:
        ones = [S[i] for i in range(N) if t[i] =='1']
        zeros = [S[i] for i in range(N) if t[i] =='0']
        X = set.intersection(*ones)
        X.difference_update(*zeros)
        weights[''.join(t)] = len(X)

        #save the intersected results
        if options.saveoverlaps:
            if len(X) >= options.overlapthresh:
                file_name = ''
                name_itr = 0
                for name in t:
                    if name == '1':
                        file_name += '_'+label_names[name_itr]
                    name_itr +=1
                file_name = ''.join(t)+file_name
                hlp.create_dir(output+'/sets')
                inter_file = open(output+'/sets/'+file_name+'.txt', 'w')
                inter_file.writelines('\n'.join(list(X)))
                inter_file.close()

    return(weights)

def create_r_script(labels, names, options):
    """
    It creates Rscript for UpSetR plot for the genomic regions.

    """
    #temp_f = tempfile.NamedTemporaryFile(delete=False)
    #temp_f = open(tempfile.mktemp(), "w")
    script_file =  options.output+'/'+str(options.project)+'_'+options.command+'.R'
    temp_f = open(script_file, 'w')
    output_name = options.output+'/'+str(options.project)+'_'+options.command+'.'+options.figtype

    temp_f.write('#!/usr/bin/env Rscript'+"\n")
    temp_f.write('if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))\n')
    temp_f.write('library("UpSetR")\n')
    if options.figtype == 'ps':
        temp_f.write('if (suppressMessages(!require("Cairo"))) suppressMessages(install.packages("Cairo", repos="http://cran.us.r-project.org"))\n')
        temp_f.write('library("Cairo")\n')
    
    if options.figtype == 'pdf' or options.figtype == 'svg':
        temp_f.write(options.figtype+'("'+output_name+'", width='+str(options.figsize[0])+', height='+str(options.figsize[1])+', onefile=FALSE, useDingbats=FALSE)'+'\n')
    
    elif options.figtype == 'ps':
        temp_f.write('cairo_ps("'+output_name+'", width='+str(options.figsize[0])+', height='+str(options.figsize[1])+')'+'\n')
    else:
        temp_f.write(options.figtype+'("'+output_name+'", width='+str(options.dpi*options.figsize[0])+', height='+str(options.dpi*options.figsize[1])+', res='+str(options.dpi)+')\n')
     
    temp_f.write("expressionInput <- c(")

    last = 1

    shiny = ""

    for key, value in labels.items(): #iteritems in python 2.7
        i = 0
        first = 1
        for x in key:
            if i == 0:
                temp_f.write("'")
      
            if x == '1':
                if first == 1:
                    temp_f.write(str(names[i]))
                    shiny += str(names[i])
                    first = 0
                else:
                    temp_f.write('&'+str(names[i]))
                    shiny += '&'+str(names[i])

            if i == len(key)-1:
                if last == len(labels):
                    temp_f.write("'="+str(value))
                    shiny += "="+str(value)

                else:
                    temp_f.write("'="+str(value)+',')
                    shiny += "="+str(value)+','
            i += 1
        last +=1
    temp_f.write(")\n")

    #options.shiny = True
    #If shiny output
    if not options.showshiny:

        shiny_import =  options.output+'/'+str(options.project)+'_'+options.command+'_combinations.txt'
        shiny_file = open(shiny_import, 'w')
        shiny_file.write("You can go to Intervene Shiny App https://asntech.shinyapps.io/Intervene-app/ and copy/paste the following intersection data to get more interactive figures.\n\n")
        shiny_file.write(shiny)
        shiny_file.close()
    
    else:
        print(shiny)

    if options.showsize:
        options.showsize = 'yes'

    #if options.ninter == 0:
    #    options.ninter = "NA"

    if not options.showzero:
        options.showzero = 'NULL'
    else:
        options.showzero = "'on'"

    temp_f.write('upset(fromExpression(expressionInput), nsets='+str(len(names))+', nintersects='+str(options.ninter)+', show.numbers="'+str(options.showsize)+'", main.bar.color="'+options.mbcolor+'", sets.bar.color="'+options.sbcolor+'", empty.intersections='+str(options.showzero)+', order.by = "'+options.order+'", number.angles = 0, mainbar.y.label ="'+options.mblabel+'", sets.x.label ="'+options.sxlabel+'")\n')
    temp_f.write('invisible(dev.off())\n')

    #print temp_f.read()
    #print temp_f.name
    #cmd = 'intervene_upset_plot.R %s %s %s' % ('genomic',5,temp_f.name)
    cmd = temp_f.name
    temp_f.close()

    if not options.scriptonly:
        os.system('chmod +x '+cmd)
        os.system(cmd)
        print('\nYou are done! Please check your results @ '+options.output+'. \nThank you for using Intervene!\n')
        sys.exit(0)
    else:
        print('\nYou are done! Please check your UpSet plot script and Shiny App input @ '+options.output+'. \nThank you for using Intervene!\n')
        sys.exit(0)

        
def draw_genomic(labels, names, output, fig_type):
    #temp_f = tempfile.NamedTemporaryFile(delete=False)
    temp_f = open(tempfile.mktemp(), "w")
    
    temp_f.write("expressionInput <- c(")
    last = 1
    for key, value in labels.items():
        i = 0
        first = 1
        for x in key:
            if i == 0:
                temp_f.write("'")      
            if x == '1':
                if first == 1:
                    temp_f.write(str(names[i]))
                    first = 0
                else:
                    temp_f.write('&'+str(names[i]))

            if i == len(key)-1:
                if last == len(labels):
                    temp_f.write("'="+str(value))
                else:
                    temp_f.write("'="+str(value)+',')
            i += 1
        last +=1
        #temp_f.write("'="+str(value)+',')
    temp_f.write(")\n")
    #print temp_f.read()
    #print temp_f.name
    temp_f.close()
    cmd = 'upset_plot_intervene.R %s %s %s %s %s ' % ('genomic',len(names),temp_f.name, output, fig_type)
    os.system(cmd)
    sys.exit(0)

def one_vs_rest_intersection(beds, peaks, output, **kwoptions):
    '''
    Compares a set of peaks with several other peaks sets.

    '''
    names = []
    matrix_file = output+'/One_vs_all_peak_set_matrix.txt'
    f = open(matrix_file, 'w')
    
    f.write('peak_id')
    #f.write('peak_id\tchrom\tstart\tend')

    for bed in beds:
        names.append(get_name(bed))
        f.write('\t' + str(get_name(bed)))
    #main_int.append(names)

    peaks = BedTool(peaks[0])
    f.write('\n')
    for i in peaks:
        #region_int = []
        peak_id = str(i.chrom)+"_"+str(i.start)+"_"+str(i.end)
        f.write(peak_id)
        
        for bed in beds:
            b = BedTool(bed)
            int_count = BedTool(str(i), from_string=True).intersect(b).count()
            if (int_count > 0):
                #region_int.append("1")
                f.write('\t' + str(1))
            else:
                #region_int.append("0")
                f.write('\t' + str(0))
        f.write('\n')
        #main_int.append(region_int)
        #matrix[peak_id] = region_int
    f.close()

    return matrix_file
    