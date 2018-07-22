#!/usr/bin/python
import csv, sys, os, errno, math, random, operator
from intervaltree import Interval, IntervalTree
from random import shuffle
import numpy as np
from copy import deepcopy

###############################################################################
#
#   CHANGE THIS
#
################################################################################

location_label = "location"         #IS location column
chrom_label = "landmark"            #chromosome column

BASE_PATH = ""                      #path to files

paths = {
    "SuperEnh":     [BASE_PATH + "SE/Jurkat.bed",
                     BASE_PATH + "SE/MNVa.bed"],
    "SE_invivo":    [BASE_PATH + "SE/MNVa.bed"],
    "SE_invitro":   [BASE_PATH + "SE/Jurkat.bed"],
    "ISGene":       [BASE_PATH + "gene.is.in.csv"],
    "IS":           [BASE_PATH + "is.in_short.csv"],    #file with IS to map
    "gene":         [BASE_PATH + "gene_list_out.csv"],
    "nongene":      [BASE_PATH + "nongenic.csv"],
    "intex":        [BASE_PATH + "exon.csv"],
    "cds":          [BASE_PATH + "cds.csv"]
}

###############################################################################


###############################################################################
#
#   map IS into genes
#
################################################################################

def map_IS(file_IS):
    new_file_IS = list(file_IS)
    is_len = len(file_IS)
    with open(paths["gene"][0], 'rt') as csvfile:
        genes = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        genes = sorted(genes, key=operator.itemgetter("landmark", "gene_start"))
        genes = list(genes)
        #create interval tree
        hash_genes = {}
        tree = {}
        i = 0
        old_chr = ''
        for gene in genes:
            chrom = gene['landmark']
            if chrom != old_chr:
                hash_genes[chrom] = {}
                hash_genes[chrom]['start'] = i
                tree[chrom] = IntervalTree()
                if old_chr != '':
                    hash_genes[old_chr]['end'] = i - 1
                old_chr = chrom
            i += 1
        hash_genes[old_chr]['end'] = i - 1
        for key in hash_genes:
            for gene_id in range(hash_genes[key]['start'], hash_genes[key]['end']+1):
                gene_start = int(genes[gene_id]['gene_start'])
                gene_end = int(genes[gene_id]['gene_end'])
                gene_chr = genes[gene_id]['landmark']
                tree[gene_chr][gene_start:gene_end] = gene_id
        jk = 0
        for issite in file_IS:
            jk += 1
            print("{}/{}".format(jk, is_len))
            #genic IS
            is_genic = False
            is_loc = int(issite[location_label])
            is_chrom = issite[chrom_label]
            point = tree[is_chrom][is_loc]
            for p in point:
                gene_id = p.data
                gene_name = genes[gene_id]['name']
                gene_start = int(genes[gene_id]['gene_start'])
                gene_end = int(genes[gene_id]['gene_end'])
                
                if is_loc in range(gene_start + 1, gene_end):
                    if is_genic == False:
                        tmp_row = issite
                    else:
                        tmp_row = issite.copy()
                    gene_len = gene_end - gene_start
                    gene_frId = is_chrom + ':' + str(gene_start) + ':' + gene_name
                    gene_orient = genes[gene_id]['gene_orientation']
                    gene_introns = genes[gene_id]['intron_count']
                    
                    if gene_orient == 'F':
                        gene_front = is_loc - gene_start
                        per_front = math.ceil(round((gene_front / gene_len)*100, 8))
                        per_front2 = gene_front/gene_len
                    else:
                        gene_front = gene_end - is_loc 
                        per_front = math.ceil(round((gene_front / gene_len)*100, 8))
                        per_front2 = gene_front/gene_len

                    tmp_row.update({'ncbi_gene_id':genes[gene_id]['ncbi_gene_id']})
                    tmp_row.update({'gene':gene_name})
                    tmp_row.update({'ncbi_gene_length':gene_len})
                    tmp_row.update({'ncbi_gene_start':gene_start})
                    tmp_row.update({'ncbi_gene_end':gene_end})
                    tmp_row.update({'ncbi_gene_orientation':gene_orient})
                    tmp_row.update({'GeneFragmentID':gene_frId})
                    tmp_row.update({'TSS_distance':gene_front})  
                    tmp_row.update({'TES_distance':gene_len - gene_front})  
                    tmp_row.update({'per_front':per_front})  
                    tmp_row.update({'per_front2':per_front2}) 
                    tmp_row.update({'intron_count':gene_introns})
                    if is_genic == True:
                        new_file_IS.append(tmp_row)
                    is_genic = True
            #integenic IS
            if is_genic == False:
                issite.update({'ncbi_gene_id':'NA'})
                issite.update({'gene':'NA'})
                issite.update({'ncbi_gene_length':'NA'})
                issite.update({'ncbi_gene_start':'NA'})
                issite.update({'ncbi_gene_end':'NA'})
                issite.update({'ncbi_gene_orientation':'NA'})
                issite.update({'GeneFragmentID':'NA'})
                issite.update({'TSS_distance':'NA'})  
                issite.update({'TES_distance':'NA'})  
                issite.update({'per_front':'NA'})  
                issite.update({'per_front2':'NA'}) 
                issite.update({'intron_count':'NA'})
    return new_file_IS      
    

###############################################################################
#
#   calculate   - if IS is inside intron/exon
#
################################################################################

def intron(file_IS):

    with open(paths["intex"][0], 'rt') as csvfile:
        genes = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        genes = sorted(genes, key=operator.itemgetter("chr", "gene"))
        genes = list(genes)
        
        #create hash table
        hash_genes = {}
        i = 0
        for gene in genes:
            if gene['type'] == 'gene':
                gene_name = gene['gene'].upper()
                gene_chrom = gene['chr']
                gene_start = gene['start']
                label = gene_chrom + ':' + gene_start + ':' + gene_name
                hash_genes[label] = i
            i += 1
        
        il = len(file_IS)
        ik = 0
        for IS in file_IS:
            ik += 1
            print("{}/{}".format(ik, il))
            is_chr = IS['landmark']
            is_geneName  = IS['gene'].upper()
            is_geneId  = IS['ncbi_gene_id']
            is_loc = int(IS['location'])
            is_gene_start = IS['ncbi_gene_start']
            is_ex = prev_ex = prev_gene = ""
            
            if is_geneId == "NA":
                IS.update({"exon/intron":""})
                IS.update({"ex_up_len":""})
                IS.update({"ex_up_type":""})
                IS.update({"ex_down_len":""})
                IS.update({"ex_down_type":""})
                IS.update({"type_len":""})
            else:

                label = is_chr + ':' + str(is_gene_start) + ':' + is_geneName
                list_id = int(hash_genes[label])
                gene_start = int(genes[list_id]['start'])
                gene_end = int(genes[list_id]['end'])
                list_id += 1
                while genes[list_id]['type'] != 'gene':
                    gene = genes[list_id]
                    gene_type = gene['type']
                    exon_gene = gene['gene'].upper()
                    #exon/intron
                    if gene_type == 'exon':
                        exon_start = int(gene['start'])
                        exon_end = int(gene['end'])
                        #print(exon_start,is_loc,exon_end)

                        #exon
                        if is_loc in range(exon_start+1, exon_end):
                            if is_ex == "":
                                is_ex = "exon"
                                is_ex_len = exon_end - exon_start
                                #distance to splice site/whether SS is upstream or downstream/whether it is donor or acceptor
                                #only one exon
                                if exon_start == gene_start and exon_end == gene_end:
                                    ex_up_len = 'NA'
                                    ex_up_type = ""
                                    ex_down_len = 'NA'
                                    ex_down_type = ""
                                elif gene['strand'] == "+":
                                    len1 = is_loc - exon_start
                                    len2 = exon_end - is_loc
                                    #first exon
                                    if exon_start == gene_start:
                                        ex_up_len = 'NA'
                                        ex_up_type = ""
                                        ex_down_len = len2
                                        ex_down_type = "SD"
                                    #last exon
                                    elif exon_end == gene_end:
                                        ex_up_len = len1
                                        ex_up_type = "SA"
                                        ex_down_len = 'NA'
                                        ex_down_type = ""
                                    else:
                                        ex_up_len = len1
                                        ex_up_type = "SA"
                                        ex_down_len = len2
                                        ex_down_type = "SD"
                                else:
                                    len1 = is_loc - exon_start
                                    len2 = exon_end - is_loc
                                    #last exon
                                    if exon_start == gene_start:
                                        ex_up_len = len2
                                        ex_up_type = "SA"
                                        ex_down_len = 'NA'
                                        ex_down_type = ""
                                    #first exon
                                    elif exon_end == gene_end:
                                        ex_up_len = 'NA'
                                        ex_up_type = ""
                                        ex_down_len = len1
                                        ex_down_type = "SD"
                                    else:
                                        ex_up_len = len2
                                        ex_up_type = "SA"
                                        ex_down_len = len1
                                        ex_down_type = "SD"
                                break
                            else:
                                string = "is {} inside two exons/introns {}".format(is_loc, is_ex)
                                print(string) 
                                sys.exit(1)
                        #intron
                        if exon_start != gene_start and is_loc in range(prev_end, exon_start+1):
                            if is_ex == "":
                                is_ex = "intron"
                                len1 = is_loc - prev_end
                                len2 = exon_start - is_loc
                                is_ex_len = exon_start - prev_end
                                #distance to splice site/whether SS is upstream or downstream/whether it is donor or acceptor
                                if gene['strand'] == "+":
                                    ex_up_len = len1
                                    ex_up_type = "SD"
                                    ex_down_len = len2
                                    ex_down_type = "SA"
                                else:
                                    ex_up_len = len2
                                    ex_up_type = "SD"
                                    ex_down_len = len1
                                    ex_down_type = "SA"
                                break
                            else:
                                string = "is {} inside two exons/introns {}".format(is_loc, is_ex)
                                print(string) 
                        prev_end = exon_end
                    list_id += 1
                IS.update({"exon/intron":is_ex})
                IS.update({"ex_up_len":ex_up_len})
                IS.update({"ex_up_type":ex_up_type})
                IS.update({"ex_down_len":ex_down_len})
                IS.update({"ex_down_type":ex_down_type})
                IS.update({"type_len":is_ex_len})

    return file_IS


###############################################################################
#
#   calculate   - if IS is inside 5'/3'UTR
#
################################################################################

def cds(file_IS):
    with open(paths["cds"][0], 'rt') as csvfile:
        cdss = csv.DictReader(csvfile, delimiter=' ', quotechar='|')
        cdss = list(cdss)
        hash_genes = {}
        i = 0
        for cds in cdss:
            gene_name = cds['gene_name'].upper()
            gene_chrom = cds['chrom']
            gene_start = cds['transcription_start']
            label = gene_chrom + ':' + gene_start + ':' + gene_name
            hash_genes[label.upper()] = i
            i += 1
        il = len(file_IS)
        ik = 0
        for IS in file_IS:
            ik += 1
            if (ik % 1000) == 0:
                print("{}/{}".format(ik, il))
            is_chr = IS['landmark']
            is_geneName  = IS['gene'].upper()
            is_geneId  = IS['ncbi_gene_id']
            is_loc = int(IS['location'])
            type2_len = -1
            if is_geneId != "NA":
                is_gene = False
                try:
                    label = hash_genes[IS['GeneFragmentID'].upper()]
                except:
                    label = -1
                    
                if label == -1:
                    IS.update({"type":""})
                    IS.update({"type2":IS['exon/intron']})
                    IS.update({"StartCodon_dist":""})
                    IS.update({"StartCodon_pos":""})
                    IS.update({"StopCodon_dist":""})
                    IS.update({"StopCodon_pos":""})
                    IS.update({"gene_strand":""})
                    IS.update({"type2_len":IS['type_len']})
                    IS.update({"translation_start":""})
                    IS.update({"translation_end":""})
                else:
                    cds = cdss[label]
                    cds_geneId = cds['gene_id']
                    cds_chr = cds['chrom']
                    transcription_start = int(cds['transcription_start'])
                    transcription_end = int(cds['transcription_end'])
                    if cds_chr == is_chr and cds_geneId == is_geneId and is_loc in range(transcription_start+1, transcription_end):
                        if is_gene:
                            print("{} inside second gene {} {}".format(is_loc, cds_geneId, cds['gene_name']))
                            sys.exit(0)
                        is_gene = True
                        if cds['strand'] == '+':
                            cds_strand = '+'
                            if cds['start'] == 'NA':
                                translation_start = -1
                            else:
                                translation_start = int(cds['start']) + int(cds['start_phase'])
                            if cds['end'] == 'NA':
                                translation_end = -1
                            else:
                                translation_end = int(cds['end'])
                            trans_start = translation_start
                            trans_end = translation_end
                            #5'UTR/3'UTR
                            if translation_start != -1 and is_loc in range(transcription_start+1, translation_start):
                                cds_type = "5'UTR"
                                type2_len = translation_start - transcription_start
                            elif translation_end != -1 and is_loc in range(translation_end+1, transcription_end):
                                cds_type = "3'UTR"
                                type2_len = transcription_end - translation_end
                            else:
                                cds_type = ""
                            #distance to start codon
                            if translation_start == -1:
                                start_dist = 'NA'
                                start_pos = 'NA'
                            elif is_loc < translation_start:
                                start_dist = translation_start - is_loc
                                start_pos = "DOWNSTREAM"
                            elif is_loc >= translation_start:
                                start_dist = is_loc - translation_start
                                start_pos = "UPSTREAM"
                            #distance to stop codon
                            if translation_end == -1:
                                stop_dist = 'NA'
                                stop_pos = 'NA'
                            elif is_loc <= translation_end:
                                stop_dist = translation_end - is_loc
                                stop_pos = "DOWNSTREAM"
                            elif is_loc > translation_end:
                                stop_dist = is_loc - translation_end
                                stop_pos = "UPSTREAM"
                        else:
                            cds_strand = '-'
                            if cds['start'] == 'NA':
                                translation_end = -1
                            else:
                                translation_end = int(cds['start'])
                            if cds['end'] == 'NA':
                                translation_start = -1
                            else:
                                translation_start = int(cds['end']) - int(cds['end_phase'])
                            transcription_start = int(cds['transcription_end'])
                            transcription_end = int(cds['transcription_start'])
                            trans_start = translation_end
                            trans_end = translation_start
                            #5'UTR/3'UTR
                            if translation_start != -1 and is_loc in range(translation_start+1, transcription_start):
                                cds_type = "5'UTR"
                                type2_len = transcription_start - translation_start
                            elif translation_end != -1 and is_loc in range(transcription_end+1, translation_end):
                                cds_type = "3'UTR"
                                type2_len = translation_end - transcription_end
                            else:
                                cds_type = ""
                            #distance to start codon
                            if translation_start == -1:
                                start_dist = 'NA'
                                start_pos = 'NA'
                            elif is_loc <= translation_start:
                                start_dist = translation_start - is_loc
                                start_pos = "UPSTREAM"
                            elif is_loc > translation_start:
                                start_dist = is_loc - translation_start
                                start_pos = "DOWNSTREAM"
                            #distance to stop codon
                            if translation_end == -1:
                                stop_dist = 'NA'
                                stop_pos = 'NA'
                            elif is_loc < translation_end:
                                stop_dist = translation_end - is_loc
                                stop_pos = "UPSTREAM"
                            elif is_loc >= translation_end:
                                stop_dist = is_loc - translation_end
                                stop_pos = "DOWNSTREAM"
                    
                    if type2_len == -1:
                        type2_len = IS['type_len']
                    if cds_type == "":
                        cds_type2 = IS['exon/intron']
                    else:
                        cds_type2 = cds_type
                    trans_start = 'NA' if trans_start == -1 else trans_start
                    trans_end = 'NA' if trans_end == -1 else trans_end
                    IS.update({"type":cds_type})
                    IS.update({"type2":cds_type2})
                    IS.update({"StartCodon_dist":start_dist})
                    IS.update({"StartCodon_pos":start_pos})
                    IS.update({"StopCodon_dist":stop_dist})
                    IS.update({"StopCodon_pos":stop_pos})
                    IS.update({"gene_strand":cds_strand})
                    IS.update({"type2_len":type2_len})
                    IS.update({"translation_start":trans_start})
                    IS.update({"translation_end":trans_end})
            else:
                IS.update({"type":""})
                IS.update({"type2":""})
                IS.update({"StartCodon_dist":""})
                IS.update({"StartCodon_pos":""})
                IS.update({"StopCodon_dist":""})
                IS.update({"StopCodon_pos":""})
                IS.update({"gene_strand":""})
                IS.update({"type2_len":IS['type_len']})
                IS.update({"translation_start":''})
                IS.update({"translation_end":''})
    return file_IS

###############################################################################
#
#   calculate   - if intron/exon is untranslated or translated
#
################################################################################   

def intron2(file_IS):

    with open(paths["intex"][0], 'rt') as csvfile:
        genes = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        genes = sorted(genes, key=operator.itemgetter("chr", "gene"))
        genes = list(genes)
        
        #create hash table
        hash_genes = {}
        i = 0
        for gene in genes:
            if gene['type'] == 'gene':
                gene_name = gene['gene'].upper()
                gene_chrom = gene['chr']
                gene_start = gene['start']
                label = gene_chrom + ':' + gene_start + ':' + gene_name
                hash_genes[label] = i
            i += 1
        
        il = len(file_IS)
        ik = 0
        for IS in file_IS:
            ik += 1
            print("{}/{}".format(ik, il))
            is_chr = IS['landmark']
            is_geneName  = IS['gene'].upper()
            is_geneId  = IS['ncbi_gene_id']
            is_loc = int(IS['location'])
            is_gene_start = IS['ncbi_gene_start']
            is_gene_orient = IS['ncbi_gene_orientation']
            is_ex = prev_ex = prev_gene = ""
            IS.update({"type3":""})
            IS.update({"type3_len":""})
            if is_geneId == "NA":
                IS.update({"type3":"intergenic"})
                IS.update({"type3_len":IS['type2_len']})
            else:
                label = is_chr + ':' + str(is_gene_start) + ':' + is_geneName
                list_id = int(hash_genes[label])
                gene_start = int(genes[list_id]['start'])
                gene_end = int(genes[list_id]['end'])
                if IS['translation_start'] == '' and IS['translation_end'] == '':
                    IS.update({"type3":IS['type2']+'T'})
                    IS.update({"type3_len":IS['type2_len']})
                    continue
                if IS['translation_start'] =='NA':
                    translation_start = gene_start - 1
                else:
                    translation_start = int(IS['translation_start'])
                if IS['translation_end'] =='NA':
                    translation_end = gene_end + 1
                else:
                    translation_end = int(IS['translation_end'])
                if translation_start > translation_end:
                    tmp = translation_start
                    translation_start = translation_end
                    translation_end = tmp
                   
                list_id += 1
                while genes[list_id]['type'] != 'gene':
                    
                    gene = genes[list_id]
                    gene_type = gene['type']
                    exon_gene = gene['gene'].upper()
                    #exon/intron
                    if gene_type == 'exon':
                        exon_start = int(gene['start'])
                        exon_end = int(gene['end'])
                        
                        #exon
                        if is_loc in range(exon_start+1, exon_end):
                            if is_ex == "":
                                #print('gene_start: {}, translation_start: {}, exon_start: {}, is_loc: {}, exon_end: {}, translation_end: {}, gene_end: {}'.format(gene_start, translation_start, exon_start, is_loc, exon_end, translation_end, gene_end))
                                if is_loc in range(gene_start+1, translation_start):
                                    type3 = 'exonU'
                                    type3_len = min(translation_start, exon_end) - max(gene_start, exon_start)
                                elif is_loc in range(translation_start, translation_end+1):
                                    type3 = 'exonT'
                                    type3_len = min(translation_end, exon_end) - max(translation_start, exon_start)
                                elif is_loc in range(translation_end+1, gene_end):
                                    type3 = 'exonU'
                                    type3_len = min(gene_end, exon_end) - max(translation_end, exon_start)
                                else:
                                    sys.exit(1)
                                break
                            else:
                                string = "is {} inside two exons/introns {}".format(is_loc, is_ex)
                                print(string) 
                                sys.exit(1)
                        #intron
                        if exon_start != gene_start and is_loc in range(prev_end, exon_start+1):
                            if is_ex == "":
                                is_ex = "intron"
                                #print('gene_start: {}, translation_start: {}, prev_end: {}, is_loc: {}, exon_start: {}, translation_end: {}, gene_end: {}'.format(gene_start, translation_start, prev_end, is_loc, exon_start, translation_end, gene_end))
                                if is_loc in range(gene_start+1, translation_start):
                                    type3 = 'intronU'
                                    type3_len = min(translation_start, exon_start) - max(gene_start, prev_end)
                                elif is_loc in range(translation_start, translation_end+1):
                                    type3 = 'intronT'
                                    type3_len = min(translation_end, exon_start) - max(translation_start, prev_end)
                                elif is_loc in range(translation_end+1, gene_end):
                                    type3 = 'intronU'
                                    type3_len = min(gene_end, exon_start) - max(translation_end, prev_end)
                                    
                                break
                            else:
                                string = "is {} inside two exons/introns {}".format(is_loc, is_ex)
                                print(string) 
                                sys.exit(1)
                        prev_end = exon_end
                    list_id += 1
                IS.update({"type3":type3})
                IS.update({"type3_len":type3_len})

    return file_IS


###############################################################################
#
#   main
#
################################################################################        
                

with open(paths["IS"][0], 'rt') as csvfile2:
    #open file with IS
    intSiteGene = csv.DictReader(csvfile2, delimiter=',', quotechar='"')
    intSiteGene = list(intSiteGene)
    in_file = list(intSiteGene)
    #map into genes
    in_file = map_IS(in_file)
    #intron exon
    in_file = intron(in_file)
    #UTR
    in_file = cds(in_file)
    #tranlated/untranslated exon/intron
    in_file = intron2(in_file)
    #save file
    with open("gene.is.in.csv", 'w') as output:
        writer = csv.writer(output, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(in_file[0])
        for row in in_file:
            writer.writerow(row.values())
