#!/usr/bin/env python3
#v0.1.5
import os
import sys
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import pdb
import re
from operator import itemgetter, attrgetter

def path_asfileslist(path):
    if os.path.isdir(path):
        files_list = []
        for (dirpath, dirnames, filenames) in os.walk(path):
            files_list.extend(filenames)
        return files_list
    else:
        raise argparse.ArgumentTypeError(f'readable_dir:{path} is not a valid path')

def string_asfileslist(string):
    files_list = []
    for path in string.split(','):
        if os.path.isdir(path):
            files_list.append(path)
        else:
            raise argparse.ArgumentTypeError(f'readable_dir:{path} is not a valid path')
    return files_list

def table_aslist(table, separator):
    table_list = []
    table_temp = [ x.split(separator) for x in open(table, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n') ]
    for i in range(len(table_temp[0])):
        table_list.append([])
    for i in range(len(table_temp)):
        for j in range(len(table_temp[i])):
            table_list[j].append(table_temp[i][j])
    return table_list

def table_asdf(table, separator):
    table_temp = [ [ y.strip('"') for y in x.split(separator)[1:] ] for x in open(table, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[1:] ]
    row_names = [ x.split(separator)[0] for x in open(table, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[1:] ]
    column_names = open(table, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[0].split(separator)[1:]
    table_df = pd.DataFrame(table_temp, columns = column_names, index=row_names)
    return table_df

def filter_table(table, filter_list):
    outtable_index = []
    outtable = []
    for i in range(len(table[0])):
        if any(f for f in filter_list if f in table[0][i]):
            outtable_index.append(i)
    for i in range(len(table)):
        outtable.append([])
    for i in range(len(table)):
       for j in range(len(outtable_index)):
           outtable[i].append(table[i][outtable_index[j]])
    return outtable

def filter_dfrownames(df, filter_list, mode):
    subsample_index = []
    for i in range(len(list(df.index))):
        if len([ x for x in filter_list if x in df.index[i]]) > 0:
            subsample_index.append(list(df.index)[i])
    if mode == 'keep':
        return df.filter(items = subsample_index, axis=0)
    elif mode == 'drop':
        return df.drop(list(df.filter(items = subsample_index, axis=0).index.values))

def gff_asposlist(gff, namefield, feature_list):
    gff_poslist = []
    gff_pnname = []
    gff_ref = []
    gff_strand = []
    gff_position = []
    for item in [ x for x in open(gff, 'r').read().replace("\r\n", "\n").rstrip('\n').split('\n') if x[0] != '#' ]:
        if item.split('\t')[2] in feature_list:
            pn = item.split('\t')
            pnname_index = [ x.split('=')[0] for x in pn[8].split(';') ].index(namefield)
            pnname = [ x.split('=')[1] for x in pn[8].split(';') ][pnname_index]

            if pnname not in gff_pnname:
                gff_pnname.append(pnname)
                gff_ref.append(pn[0])
                gff_strand.append(pn[6])
                gff_position.append([])
            
            if gff_strand[gff_pnname.index(pnname)] == '+':
                gff_position[gff_pnname.index(pnname)] += [ x for x in range(int(pn[3]),int(pn[4])+1) ]
            elif gff_strand[gff_pnname.index(pnname)] == '-':
                gff_position[gff_pnname.index(pnname)] += [ -x for x in range(int(pn[4]),int(pn[3])-1,-1) ]
    for i in range(len(gff_pnname)):
        gff_poslist.append([gff_pnname[i], gff_position[i], [x+1 for x in range(int(len(gff_position[i])/3)) if len(gff_position[i])%3 == 0]])
    return gff_poslist

def seek_dp(value, separator, threshold):
    value_list = [ int(x) for x in str(value).split(separator) if x.isdigit() ]
    if len(value_list) == 0:
        return "NA"
    elif max(value_list) >= threshold:
        return "DP"
    else:
        return "NO"

def seek_varcount(value, separator, threshold):
    value_list = [ int(x) for x in str(value).split(separator) if x.isdigit() ]
    if np.isnan(value):
        return "NA"
    elif int(value) >= threshold:
        return "VARCOUNT" + str(int(value))
    else:
        return "NO"
        
def seek_poi(var, value, poi):

    if len(value) == 0 or not value[0].isdigit() or poi[0][0] == '':
        return ""

    poi_var = poi[0][:]
    poi_comment = poi[1][:]
    poi_comment.append('frameshift')
    poi_comment_set = list(set(poi_comment))
    var_found = [ var[i] for i in range(len(var)) if int(value[i]) > 0 ]
    var_found_maj = [ var[i] for i in range(len(var)) if int(value[i]) >= 50 ]
    value_found = [ value[i] for i in range(len(var)) if int(value[i]) > 0 ]
    value_found_maj = [ value[i] for i in range(len(var)) if int(value[i]) >= 50 ]

    poi_found = []
    for i in range(len(poi_var)):
        for j in range(len(var_found)):
            if poi_var[i] in var_found[j]:
                poi_found.append([poi_comment[i], var_found[j], value_found[j]])

    for i in range(len(var_found_maj)):
        if (var_found_maj[i].split('|')[6] == 'DEL' or var_found_maj[i].split('|')[6] == 'INS' or var_found_maj[i].split('|')[6] == 'COMPLEX') and var_found_maj[i].split('|')[2] == 'CODING':
            if abs(len(var_found_maj[i].split('|')[7].split('>')[0]) - len(var_found_maj[i].split('|')[7].split('>')[1])) % 3 > 0 and ( var_found_maj[i].split('|')[4] == 'HA' or var_found_maj[i].split('|')[4] == 'NA' ):
                poi_found.append(['frameshift', var_found_maj[i], value_found_maj[i]])

    poi_formatted = []
    for i in range(len(poi_found)):
        type_temp = poi_found[i][1].split('|')[6]
        if type_temp == 'DEL':
            poi_formatted.append([poi_found[i][0], poi_found[i][1].split('|')[1] + ':' + poi_found[i][1].split('|')[7].split('>')[0][1:] + str(int(poi_found[i][1].split('|')[3])+1) + '-' + '_' + poi_found[i][2] + 'perc']) 
        elif type_temp == 'INS':
            poi_formatted.append([poi_found[i][0], poi_found[i][1].split('|')[1] + ':' + 'INS' + str(int(poi_found[i][1].split('|')[3])+1) + poi_found[i][1].split('|')[7].split('>')[1][1:] + '_' + poi_found[i][2] + 'perc']) 
        elif type_temp == 'COMPLEX':
            poi_formatted.append([poi_found[i][0], poi_found[i][1].split('|')[1] + ':' + 'COMPLEX' + poi_found[i][1].split('|')[7].split('>')[0] + str(int(poi_found[i][1].split('|')[3])) + poi_found[i][1].split('|')[7].split('>')[1] + '_' + poi_found[i][2] + 'perc']) 
        elif type_temp == 'MNP':
            poi_formatted.append([poi_found[i][0], poi_found[i][1].split('|')[1] + ':' + 'MNP' + poi_found[i][1].split('|')[7].split('>')[0] + str(int(poi_found[i][1].split('|')[3])) + poi_found[i][1].split('|')[7].split('>')[1] + '_' + poi_found[i][2] + 'perc']) 
        else:
            poi_formatted.append([poi_found[i][0], poi_found[i][1].split('|')[4] + ':' + poi_found[i][1].split('|')[6].split('>')[0] + str(int(poi_found[i][1].split('|')[5])) + poi_found[i][1].split('|')[6].split('>')[1] + '_' + poi_found[i][2] + 'perc']) 
    
    poi_final = []
    for i in range(len(poi_comment_set)):
        poi_regroup = []
        for x in poi_formatted:
            if x[0] == poi_comment_set[i]:
                poi_regroup.append(x[1])
        if len(poi_regroup) > 0 :
            poi_final.append(poi_comment_set[i] + ' : ' + ';'.join(poi_regroup))

    return '; '.join(poi_final)
        
def gen_metrics(value, separator, threshold):
    value_list = [ int(x) for x in str(value).split(separator) if x.isdigit() ]
    if len(value_list) == 0:
        return "NA"
    elif max(value_list) >= threshold:
        return "VARCOUNT>=" + str(threshold) 
    else:
        return "NO"

def validate_column_ncov(value, criteria1, criteria2, separator, threshold, qc):
    criteria1_list = [ float(x) for x in str(criteria1).split(separator) if x != 'NA' ]
    criteria2_list = [ float(x)*100 for x in str(criteria2).split(separator) if x != 'NA' ]
    if value in ['', 'nan'] and min(criteria1_list) >= threshold and min(criteria2_list) >= threshold:
        return "ABSCLADE"
    elif len(criteria1_list) == 0 or len(criteria2_list) == 0:
        return "ININT"
    elif min(criteria1_list) < threshold or min(criteria2_list) < threshold:
        return "ININT"
    else:
        return value

def validate_column_fluabv(value, criteria1, criteria2, separator, threshold, qc):
    criteria1_list = [ float(x) for x in str(criteria1).split(separator) if x != 'NA' ]
    criteria2_list = [ float(x) for x in str(criteria2).split(separator) if x != 'NA' ]
    if qc not in ['good', 'mediocre']:
        return "ININT"
    elif value in ['', 'nan'] and min(criteria1_list) >= threshold and min(criteria2_list) >= threshold:
        return "ABSCLADE"
    elif len(criteria1_list) == 0 or len(criteria2_list) == 0:
        return "ININT"
    elif min(criteria1_list) < threshold or min(criteria2_list) < threshold:
        return "ININT"
    else:
        return value

def validate_column_hmpxv(value, criteria1, separator, threshold, qc):
    criteria1_list = [ float(x) for x in str(criteria1).split(separator) if x != 'NA' ]
    if qc not in ['good', 'mediocre']:
        return "ININT"
    elif value in ['', 'nan'] and min(criteria1_list) >= threshold:
        return "ABSCLADE"
    elif len(criteria1_list) == 0:
        return "ININT"
    elif min(criteria1_list) < threshold:
        return "ININT"
    else:
        return value

def validate_column_rsv(value, criteria1, criteria2, criteria3, separator, threshold, qc):
    criteria1_list = [ float(x) for x in str(criteria1).split(separator) if x != 'NA' ]
    criteria2_list = [ float(x)*100 for x in str(criteria2).split(separator) if x != 'NA' ]
    criteria3_list = [ float(x)*100 for x in str(criteria3).split(separator) if x != 'NA' ]
    if qc not in ['good', 'mediocre', 'bad']:
        return "ININT"
    elif value in ['', 'nan'] and min(criteria1_list) >= threshold and min(criteria2_list) >= threshold and min(criteria3_list) >= threshold:
        return "ABSCLADE"
    elif len(criteria1_list) == 0 or len(criteria2_list) == 0 or len(criteria3_list) == 0:
        return "ININT"
    elif min(criteria1_list) < threshold or min(criteria2_list) < threshold or min(criteria3_list) < threshold:
        return "ININT"
    else:
        return value

def list_position(value, valuesep, rangesep):
    if value == '' or value == 'nan':
        poslist = []
    else:
        poslist = []
        rangelist = [ x for x in value.split(rangesep) ]
        for x in rangelist:
            poslist.append(range(int(x.split(valuesep)[0]),int(x.split(valuesep)[-1])+1))
        poslist = [pos for plist in poslist for pos in plist]
    return poslist

def list_position_fluabv(value, valuesep, rangesep, insertions):
    if (value == '' or value == 'nan') and insertions == '':
        poslist = []
    else:
        insertions_list = insertions.split(',')
        if len(insertions_list) > 0:
            inspos = [ x.split(':')[0] for x in insertions_list]
        else:
            inspos = []
        if '0' in inspos:
            prefixlen = len(insertions_list[inspos.index('0')].split(':')[1])
        else:
            prefixlen = 0
        poslist = []
        rangelist = [ x for x in value.split(rangesep) if x != 'nan' ]
        for x in rangelist:
            poslist.append(range(int(x.split(valuesep)[0])+prefixlen,int(x.split(valuesep)[-1])+1+prefixlen))
        poslist = [pos for plist in poslist for pos in plist]
    return poslist

def list_intersect(containedl, containingl, mode, threshold):
    
    containedl = [ str(x) for x in containedl ]
    containingl = [ str(x) for x in containingl ]
    
    if mode == 'match':
        return ';'.join([ x for x in containedl if x in containingl ])
    elif mode == 'missing':
        return ';'.join([ x for x in containedl if x not in containingl ])
    
    
    foundpos = [ x for x in containedl if x in containingl ]
    
    nbcontained = len(containingl)-len(foundpos)
    perccontained = float(nbcontained)/len(containingl)
    if mode == 'percent':
        return round(perccontained, 3)
    if mode == 'include':
        criteria = perccontained
    elif mode == 'exclude':
        criteria = 100-perccontained
    if criteria < threshold:
        return "NO"
    else:
        return "YES"

def seek_missingpos(missingabspos, gffposlist, mode):
    missingpos = []
    if mode == 'aa':
        for i in range(len(gffposlist)):
            for x in gffposlist[i][2]:
                loop_index = gffposlist[i][2].index(x)
                if gffposlist[i][1][loop_index*3] in missingabspos or gffposlist[i][1][loop_index*3+1] in missingabspos or gffposlist[i][1][loop_index*3+2] in missingabspos:
                    missingpos.append([ gffposlist[i][0] + ':' + str(x) ])
    elif mode == 'nt':
        for i in range(len(gffposlist)):
            missingpos.append([ gffposlist[i][0] + ':' + str(gffposlist[i][1].index(x)+1) for x in gffposlist[i][1] if x in missingabspos ])
    missingpos = [pos for plist in missingpos for pos in plist]
    return missingpos

def seek_missingpos_fluabv(missingabspos, gffposlist, mode):
    missingpos = []
    if mode == 'aa':
        for i in range(len(gffposlist)):
            if gffposlist[i][0] == 'SigPep' or gffposlist[i][0] == 'HA1' or gffposlist[i][0] == 'HA2':
                for x in gffposlist[i][2]:
                    loop_index = gffposlist[i][2].index(x)
                    if gffposlist[i][1][loop_index*3] in missingabspos or gffposlist[i][1][loop_index*3+1] in missingabspos or gffposlist[i][1][loop_index*3+2] in missingabspos:
                        missingpos.append([ gffposlist[i][0] + ':' + str(x) ])
    elif mode == 'nt':
        for i in range(len(gffposlist)):
            missingpos.append([ gffposlist[i][0] + ':' + str(gffposlist[i][1].index(x)+1) for x in gffposlist[i][1] if x in missingabspos ])
    missingpos = [pos for plist in missingpos for pos in plist]
    return missingpos

def match_expectedmatrix_ncov(Gclade, clade, header, mode):
    if len([ x for x in header[4:] if len(x.split('_'))!= 3 ]) > 0:
        sys.exit("expectedmatrix seems to be malformed : " + ';'.join([ x for x in header[4:] if len(x.split('_'))!= 3 ]))
    header_clade = [ x.split('_')[0] for x in header]
    header_comment_temp = [ x.split('_')[1] if len(x.split('_'))== 3 else 'NA' for x in header]
    header_comment = [ x if x!='' else '' for x in header_comment_temp]
    header_Gclade = [ x.split('_')[-1] for x in header]
    if Gclade in header_Gclade:
        if mode == 'match':
            return Gclade
        if mode == 'index':
            return header_Gclade.index(Gclade)
        if mode == 'comment':
            return header_comment[header_Gclade.index(Gclade)]
    elif '.'.join(Gclade.split('.')[0:-1]) in header_Gclade and '.'.join(Gclade.split('.')[0:-1]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-1])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))]
    elif '.'.join(Gclade.split('.')[0:-2]) in header_Gclade and '.'.join(Gclade.split('.')[0:-2]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-2])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))]
    elif '.'.join(Gclade.split('.')[0:-3]) in header_Gclade and '.'.join(Gclade.split('.')[0:-3]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-3])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))]
    elif clade in header_clade:
        if mode == 'match':
            return clade
        if mode == 'index':
            return header_clade.index(clade)
        if mode == 'comment':
            return header_comment[header_clade.index(clade)]
    else:
        if mode == 'match':
            return 'NA'
        if mode == 'index':
            return header_clade.index('ININT')
        if mode == 'comment':
            return ''

def match_expectedmatrix_fluabv(clade, root, header, mode):
    if len([ x for x in header[4:] if len(x.split('_'))!= 3 ]) > 0:
        sys.exit("expectedmatrix seems to be malformed : " + ';'.join([ x for x in header[4:] if len(x.split('_'))!= 3 ]))
    header_root = [ x.split('_')[0] for x in header]
    header_comment_temp = [ x.split('_')[1] if len(x.split('_'))== 3 else 'NA' for x in header]
    header_comment = [ x if x!='' else '' for x in header_comment_temp]
    header_clade = [ x.split('_')[-1] for x in header]
    if clade in header_clade:
        if mode == 'match':
            return clade
        if mode == 'index':
            return header_clade.index(clade)
        if mode == 'comment':
            return header_comment[header_clade.index(clade)]
    elif '.'.join(clade.split('.')[0:-1]) in header_clade and '.'.join(clade.split('.')[0:-1]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(clade.split('.')[0:-1])
        if mode == 'index':
            return header_clade.index('.'.join(clade.split('.')[0:-1]))
        if mode == 'comment':
            return header_comment[header_clade.index('.'.join(clade.split('.')[0:-1]))]
    elif '.'.join(clade.split('.')[0:-2]) in header_clade and '.'.join(clade.split('.')[0:-2]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(clade.split('.')[0:-2])
        if mode == 'index':
            return header_clade.index('.'.join(clade.split('.')[0:-2]))
        if mode == 'comment':
            return header_comment[header_clade.index('.'.join(clade.split('.')[0:-2]))]
    elif '.'.join(clade.split('.')[0:-3]) in header_clade and '.'.join(clade.split('.')[0:-3]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(clade.split('.')[0:-3])
        if mode == 'index':
            return header_clade.index('.'.join(clade.split('.')[0:-3]))
        if mode == 'comment':
            return header_comment[header_clade.index('.'.join(clade.split('.')[0:-3]))]
    else:
        if mode == 'match':
            return 'NA'
        if mode == 'index':
            return header_root.index('ININT')
        if mode == 'comment':
            return ''

def match_expectedmatrix_hmpxv(Gclade, clade, header, mode):
    if len([ x for x in header[4:] if len(x.split('_'))!= 3 ]) > 0:
        sys.exit("expectedmatrix seems to be malformed : " + ';'.join([ x for x in header[4:] if len(x.split('_'))!= 3 ]))
    header_clade = [ x.split('_')[0] for x in header]
    header_comment_temp = [ x.split('_')[1] if len(x.split('_'))== 3 else 'NA' for x in header]
    header_comment = [ x if x!='' else '' for x in header_comment_temp]
    header_Gclade = [ x.split('_')[-1] for x in header]
    if Gclade in header_Gclade:
        if mode == 'match':
            return Gclade
        if mode == 'index':
            return header_Gclade.index(Gclade)
        if mode == 'comment':
            return header_comment[header_Gclade.index(Gclade)]
    elif '.'.join(Gclade.split('.')[0:-1]) in header_Gclade and '.'.join(Gclade.split('.')[0:-1]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-1])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))]
    elif '.'.join(Gclade.split('.')[0:-2]) in header_Gclade and '.'.join(Gclade.split('.')[0:-2]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-2])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))]
    elif '.'.join(Gclade.split('.')[0:-3]) in header_Gclade and '.'.join(Gclade.split('.')[0:-3]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-3])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))]
    elif clade in header_clade:
        if mode == 'match':
            return clade
        if mode == 'index':
            return header_clade.index(clade)
        if mode == 'comment':
            return header_comment[header_clade.index(clade)]
    else:
        if mode == 'match':
            return 'NA'
        if mode == 'index':
            return header_clade.index('ININT')
        if mode == 'comment':
            return ''

def match_expectedmatrix_rsv(Gclade, clade, header, mode):
    if len([ x for x in header[4:] if len(x.split('_'))!= 3 ]) > 0:
        sys.exit("expectedmatrix seems to be malformed : " + ';'.join([ x for x in header[4:] if len(x.split('_'))!= 3 ]))
    header_clade = [ x.split('_')[0] for x in header]
    header_comment_temp = [ x.split('_')[1] if len(x.split('_'))== 3 else 'NA' for x in header]
    header_comment = [ x if x!='' else '' for x in header_comment_temp]
    header_Gclade = [ x.split('_')[-1] for x in header]
    if Gclade in header_Gclade:
        if mode == 'match':
            return Gclade
        if mode == 'index':
            return header_Gclade.index(Gclade)
        if mode == 'comment':
            return header_comment[header_Gclade.index(Gclade)]
    elif '.'.join(Gclade.split('.')[0:-1]) in header_Gclade and '.'.join(Gclade.split('.')[0:-1]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-1])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-1]))]
    elif '.'.join(Gclade.split('.')[0:-2]) in header_Gclade and '.'.join(Gclade.split('.')[0:-2]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-2])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-2]))]
    elif '.'.join(Gclade.split('.')[0:-3]) in header_Gclade and '.'.join(Gclade.split('.')[0:-3]) not in ['','A','B','C']:
        if mode == 'match':
            return '.'.join(Gclade.split('.')[0:-3])
        if mode == 'index':
            return header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))
        if mode == 'comment':
            return header_comment[header_Gclade.index('.'.join(Gclade.split('.')[0:-3]))]
    elif clade in header_clade:
        if mode == 'match':
            return clade
        if mode == 'index':
            return header_clade.index(clade)
        if mode == 'comment':
            return header_comment[header_clade.index(clade)]
    else:
        if mode == 'match':
            return 'NA'
        if mode == 'index':
            return header_clade.index('ININT')
        if mode == 'comment':
            return ''

def seek_profile(poslist, aasub):
    profile = []
    for i in range(len(poslist)):
        if aasub[i] != '':
            profile.append(poslist[i] + aasub[i])
    return profile

def seek_nocovsub(poslist, missingpos):
    #expectedpos_nt = [ [(x.split(':')[0] + ':' + str((int(''.join(filter(str.isdigit, x.split(':')[-1]))))*3-2)),(x.split(':')[0] + ':' + str((int(''.join(filter(str.isdigit, x.split(':')[-1]))))*3-1)),(x.split(':')[0] + ':' + str((int(''.join(filter(str.isdigit, x.split(':')[-1]))))*3))] for x in poslist if x.split(':')[0] != 'Ins' ]
    expectedpos_aa = [ (x.split(':')[0] + ':' + str((int(''.join(filter(str.isdigit, x.split(':')[-1])))))) for x in poslist if x.split(':')[0] != 'Ins' ]
    nocovsub = []
    for i in range(len(expectedpos_aa)):
       if expectedpos_aa[i] in missingpos:
           nocovsub.append(expectedpos_aa[i])
    return nocovsub

def seek_nocovins(poslist, missingpos):
    expectedpos_nt = [ str((int(''.join(filter(str.isdigit, x.split(':')[-1]))))) for x in poslist if x.split(':')[0] == 'Ins' ]
    expectedpos = [ x.split(':')[1] for x in poslist if x.split(':')[0] == 'Ins' ]
    nocovins = []
    for i in range(len(expectedpos)):
       if int(expectedpos_nt[i]) in missingpos:
           nocovins.append(expectedpos[i])
    return nocovins

def seek_insertion(insstring):
    inslist = [ x for x in insstring.split(',') if x not in ['', 'nan'] ]
    insertions = []
    for i in range(0, len(inslist)):
        if inslist[i].split(':')[0] != '0' and inslist[i].split(':')[1][-1] != 'N':
            insertions.append('Ins:' + inslist[i].split(':')[0] + inslist[i].split(':')[1])
            insertions.append('InsAVISBIO')
    return ','.join(insertions)

def seek_insertion_ncov(insstring):
    inslist = [ x for x in insstring.split(',') if x not in ['', 'nan'] ]
    insertions = []
    for i in range(0, len(inslist)):
        if inslist[i] == '22204:GAGCCAGAA':
            insertions.append('Ins:22204GAGCCAGAA')
        elif int(inslist[i].split(':')[0]) in list_nts:
            insertions.append('Ins:' + inslist[i].split(':')[0] + inslist[i].split(':')[1])
    return ','.join(insertions)

def seek_insertion_rsv(insstring):
    inslist = [ x for x in insstring.split(',') if x not in ['', 'nan'] ]
    insertions = []
    for i in range(0, len(inslist)):
        inslength = len((''.join(filter(lambda x: x.isdigit == False, inslist[i].split(':')[1]))))
        if inslist[i].split(':')[0] in list_nt_cds or str(int(inslist[i].split(':')[0])+inslength) in list_nt_cds:
            insertions.append('Ins:' + inslist[i].split(':')[0] + inslist[i].split(':')[1])
            insertions.append('InsAVISBIO')
    return ','.join(insertions)

def seek_issue_ncov(glims, nextclade, qc_seqcontrol, val_spikecoverage, val_atypicsub, val_atypicindel, val_nocovsub, val_nocovins, val_missingsub, nextclade_qc_privateMutations_status, summary_vcf_coinf01score, summary_vcf_coinf02match, summary_vcf_coinf02count, summary_vcf_coinf02score, val_poi, val_verif, coinf_threshold):
    result = ''
    if val_nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovsub ]
    else:
        nocovsubpos = []
    if val_nocovins != []:
        nocovsubins = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovins ]
    else:
        nocovsubins = []

    if nextclade == 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if qc_seqcontrol == 'FAILED':
            result += '_SEQFAILED'

    if nextclade != 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        ##if len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('S:')[:-1] not in [''] and x.split(':')[0] == 'S'], nocovsubpos, 'missing', 100)) >=1 or len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('Ins:')[:-1] not in [''] and x.split(':')[0] == 'Ins'], nocovsubins, 'missing', 100)) >= 1:
        ##    result += '_MISSINGSUB'
        if nextclade_qc_privateMutations_status != "good" and nextclade_qc_privateMutations_status != "mediocre":
            result += '_PM'
        if (float(summary_vcf_coinf01score) >= 1.00 and float(summary_vcf_coinf02count) >= 2 and float(summary_vcf_coinf02score) >= coinf_threshold):
            result += '_COINF-' + summary_vcf_coinf02match.replace('_', '.') + '-' + summary_vcf_coinf02score
        if val_poi != '':
            result += '_POI'
    if result != '':
        result = 'AVISBIO' + result
    return result

def seek_issue_fluabv(glims, nextclade, qc_seqcontrol, val_atypicsub, val_atypicindel, val_nocovsub, val_nocovins, val_missingsub, nextclade_qc_privateMutations_status, summary_vcf_coinf01score, summary_vcf_coinf02match, summary_vcf_coinf02count, summary_vcf_coinf02score, val_poi, val_verif, coinf_threshold):
    result = ''
    if val_nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovsub ]
    else:
        nocovsubpos = []
    if val_nocovins != []:
        nocovsubins = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovins ]
    else:
        nocovsubins = []

    if nextclade == 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if qc_seqcontrol == 'FAILED':
            result += '_SEQFAILED'
    if nextclade != 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.split(':')[-1] not in [''] ], nocovsubpos, 'missing', 100)) + len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('Ins:')[:-1] not in [''] and x.split(':')[0] == 'Ins'], nocovsubins, 'missing', 100)) >= 2:
            result += '_MISSINGSUB'
        if nextclade_qc_privateMutations_status != "good" and nextclade_qc_privateMutations_status != "mediocre":
            result += '_PM'
        if (float(summary_vcf_coinf01score) >= 1.00 and float(summary_vcf_coinf02count) >= 2 and float(summary_vcf_coinf02score) >= coinf_threshold):
            result += '_COINF-' + summary_vcf_coinf02match.replace('_', '.') + '-' + summary_vcf_coinf02score
        if val_poi != '':
            result += '_POI'
        if val_verif != 'OK':
            result += '_VERIF'
    if result != '':
        result = 'AVISBIO' + result
    return result

def seek_issue_hmpxv(glims, nextclade, qc_seqcontrol, val_atypicsub, val_atypicindel, val_nocovsub, val_nocovins, val_missingsub, nextclade_qc_privateMutations_status, summary_vcf_coinf01score, summary_vcf_coinf02match, summary_vcf_coinf02count, summary_vcf_coinf02score, val_poi, val_verif, coinf_threshold):
    result = ''
    if val_nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovsub ]
    else:
        nocovsubpos = []
    if val_nocovins != []:
        nocovsubins = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovins ]
    else:
        nocovsubins = []

    if nextclade == 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if qc_seqcontrol == 'FAILED':
            result += '_SEQFAILED'
    if nextclade != 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.split(':')[-1] not in [''] ], nocovsubpos, 'missing', 100) != '' or list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('Ins:')[:-1] not in [''] and x.split(':')[0] == 'Ins'], nocovsubins, 'missing', 100) != '':
            result += '_MISSINGSUB'
        if nextclade_qc_privateMutations_status != "good":
            result += '_PM'
        if (float(summary_vcf_coinf01score) >= 1.00 and float(summary_vcf_coinf02count) >= 2 and float(summary_vcf_coinf02score) >= coinf_threshold):
            result += '_COINF-' + summary_vcf_coinf02match.replace('_', '.') + '-' + summary_vcf_coinf02score
        if val_poi != '':
            result += '_POI'
    if result != '':
        result = 'AVISBIO' + result
    return result

def seek_issue_rsv(glims, nextclade, qc_seqcontrol, val_Fcoverage, val_Gcoverage, val_atypicsub, val_atypicindel, val_nocovsub, val_nocovins, val_missingsub, nextclade_qc_privateMutations_status, summary_vcf_coinf01score, summary_vcf_coinf02match, summary_vcf_coinf02count, summary_vcf_coinf02score, val_poi, val_verif, coinf_threshold):
    result = ''
    if val_nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovsub ]
    else:
        nocovsubpos = []
    if val_nocovins != []:
        nocovsubins = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_nocovins ]
    else:
        nocovsubins = []

    if nextclade == 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if qc_seqcontrol == 'FAILED':
            result += '_SEQFAILED'

    if nextclade != 'ININT':
        if df_validation['val_glims'].to_list().count(glims) > 1:
            result += '_DUPLICATE'
        if len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('G:')[:-1] not in [''] and x.split(':')[0] == 'G'], nocovsubpos, 'missing', 100)) >=1 or len(list_intersect([ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in val_missingsub.split(';') if x.lstrip('Ins:')[:-1] not in [''] and x.split(':')[0] == 'Ins'], nocovsubins, 'missing', 100)) >= 1:
            result += '_MISSINGSUB'
        if nextclade_qc_privateMutations_status == "":
            result += '_PM'
        if (float(summary_vcf_coinf01score) >= 1.00 and float(summary_vcf_coinf02count) >= 2 and float(summary_vcf_coinf02score) >= coinf_threshold):
            result += '_COINF-' + summary_vcf_coinf02match.replace('_', '.') + '-' + summary_vcf_coinf02score
#        if val_poi != '':
#            result += '_POI'
        if val_verif != 'OK':
            result += '_VERIF'
    if result != '':
        result = 'AVISBIO' + result
    return result

def validate_result(clade, dp, varcount, qc_seqcontrol, val_avisbio):
    result = ''
    
    if clade == 'ININT':
        if qc_seqcontrol == 'FAILED':
            return 'REPASSE SEQ FAILED'
        else:
            return 'ININT'
    else:
        result = clade
    if (dp != 'NO' and dp != 'NA') or (varcount != 'NO' and varcount != 'NA'):
        result += '_REPASSE'
    if dp != 'NO' and dp != 'NA':
        result += '_' + dp
    if varcount != 'NO' and varcount != 'NA':
        result += '_' + varcount
    if val_avisbio != '':
        result += '_' + val_avisbio
    result += ' ' + args.nextcladeversion
    return result

def validate_result_fluabv(clade, dp, varcount, qc_seqcontrol, val_avisbio):
    result = ''

    if clade == 'ININT':
        if qc_seqcontrol == 'FAILED':
            return 'REPASSE SEQ FAILED'
        else:
            return 'ININT'
    else:
        result = datatype + ' - Clade ' + clade
    if (dp != 'NO' and dp != 'NA') or (varcount != 'NO' and varcount != 'NA'):
        result += '_REPASSE'
    if dp != 'NO' and dp != 'NA':
        result += '_' + dp
    if varcount != 'NO' and varcount != 'NA':
        result += '_' + varcount
    if val_avisbio != '':
        result += '_' + val_avisbio
    result += ' ' + args.nextcladeversion
    return result

def validate_result_hmpxv(clade, Gclade, dp, varcount, qc_seqcontrol, val_avisbio):
    result = ''
    
    if clade == 'ININT':
        if qc_seqcontrol == 'FAILED':
            return 'REPASSE SEQ FAILED'
        else:
            return 'ININT'
    else:
        result = 'Clade ' + clade + ' - Gclade ' + Gclade
    if (dp != 'NO' and dp != 'NA') or (varcount != 'NO' and varcount != 'NA'):
        result += '_REPASSE'
    if dp != 'NO' and dp != 'NA':
        result += '_' + dp
    if varcount != 'NO' and varcount != 'NA':
        result += '_' + varcount
    if val_avisbio != '':
        result += '_' + val_avisbio
    result += ' ' + args.nextcladeversion
    return result

def validate_result_rsv(clade, dp, varcount, qc_seqcontrol, val_avisbio):
    result = ''
    
    if clade == 'ININT':
        if qc_seqcontrol == 'FAILED':
            return 'REPASSE SEQ FAILED'
        else:
            return 'ININT'
    else:
        result = datatype.upper() + ' - Clade ' + clade
    if (dp != 'NO' and dp != 'NA') or (varcount != 'NO' and varcount != 'NA'):
        result += '_REPASSE'
    if dp != 'NO' and dp != 'NA':
        result += '_' + dp
    if varcount != 'NO' and varcount != 'NA':
        result += '_' + varcount
    if val_avisbio != '':
        result += '_' + val_avisbio
    result += ' ' + args.nextcladeversion
    return result

def gen_commentary_ncov(classmatch, substitutions, deletions, insertions, nocovsub, nocovins, comment, atypicsub, atypicindel, val_poi):
    if classmatch == 'ININT':
        return 'DSC'
    commentary = 'Profil de la spike : '
    commentary += ';'.join(substitutions)
    if len(deletions) > 0:
         commentary += '|' + ';'.join(deletions)
    if len(insertions) > 0:
        commentary += '|' + ';'.join(insertions)
    if nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in nocovsub ]
        if len(list_intersect(nocovsubpos, list_aarbm, 'match', 100)) > 0:
            commentary += " ~ Rendu sous reserve compte tenu de(s) position(s) d'interet suivante(s) non couvertes sur la spike : " + ';'.join(nocovsub)
        else:
            commentary += " ~ Position(s) d'interet suivante(s) non couvertes sur la spike : " + ';'.join(nocovsub)
    commentary += " ~ " + comment + " (" + expectedmatrix_version[4:] + "." + expectedmatrix_version[2:4] + "." + expectedmatrix_version[0:2] + ")"
    if (atypicsub != [''] or atypicindel != ['']) and classmatch not in list_expectedmatrix_header_clade:
        commentary += " ~ Profil rare :"
        if atypicsub != ['']:
            commentary += " substitution(s) " + ";".join(atypicsub) + " presente(s) sur le RBD"
        if atypicsub != [''] and atypicindel != ['']:
            commentary += " -"
        if atypicindel != ['']:
            commentary += " indel(s) " + ";".join(atypicindel) + " present(s)"
    elif classmatch not in list_expectedmatrix_header_clade:
        commentary += ""
    if val_poi != '':
        commentary += " ~ " + val_poi
    commentary += '.'
    return commentary

def gen_commentary_fluabv(classmatch, substitutions, deletions, insertions, nocovsub, nocovins, comment, atypicsub, atypicindel, val_poi):
    if classmatch == 'ININT':
        return 'DSC'
    if comment != '':
        commentary = comment + " (Hemisphere Nord 2023)"
    else:
        commentary = "Clade distinct de ceux des souches de la composition vaccinale (Hemisphere Nord 2023)"
    #commentary = "Profil des positions d'interet de S4 : "
    #commentary += ';'.join(substitutions)
    #if len(deletions) > 0:
    #     commentary += '|' + ';'.join(deletions)
    #if len(insertions) > 0:
    #    commentary += '|' + ';'.join(insertions)
    # Ins not covered not considered
    #if nocovsub != []:
    #    nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in nocovsub ]
    #    nocovsubpos.sort()
    #    commentary += " ~ Position(s) d'interet suivante(s) non couvertes sur S4 : " + ';'.join(nocovsub)
    #if (atypicsub != [''] or atypicindel != ['']):
    #    commentary += " ~ Profil rare :"
    #    if atypicsub != ['']:
    #        commentary += " substitution(s) " + ";".join(atypicsub) + " presente(s) sur S4"
    #    if atypicsub != [''] and atypicindel != ['']:
    #        commentary += " -"
    #    if atypicindel != ['']:
    #        commentary += " indel(s) " + ";".join(atypicindel) + " present(s)"
    #elif classmatch not in list_expectedmatrix_header_clade:
    #    commentary += ""
    if val_poi != '':
        commentary += " ~ " + val_poi
    commentary += '.'
    return commentary

def gen_commentary_hmpxv(classmatch, substitutions, deletions, insertions, nocovsub, nocovins, comment, atypicsub, atypicindel, val_poi):
    if classmatch == 'ININT':
        return 'DSC'
    commentary = "Profil des positions d'interet : "
    commentary += ';'.join(substitutions)
    if len(deletions) > 0:
         commentary += '|' + ';'.join(deletions)
    if len(insertions) > 0:
        commentary += '|' + ';'.join(insertions)
    # Ins not covered not considered
    if nocovsub != []:
        nocovsubpos = [ (int(''.join(filter(str.isdigit, x.split(':')[-1])))) for x in nocovsub ]
        nocovsubpos.sort()
        commentary += " ~ Position(s) d'interet suivante(s) non couvertes : " + ';'.join(nocovsub)
    if comment != '':
        commentary += " ~ " + comment + " (" + expectedmatrix_version[4:] + "." + expectedmatrix_version[2:4] + "." + expectedmatrix_version[0:2] + ")"
    if (atypicsub != [''] or atypicindel != ['']) and atypicsub == ['TEST']:
        commentary += " ~ Profil rare :"
        if atypicsub != ['']:
            commentary += " substitution(s) " + ";".join(atypicsub) + " presente(s)"
        if atypicsub != [''] and atypicindel != ['']:
            commentary += " -"
        if atypicindel != ['']:
            commentary += " indel(s) " + ";".join(atypicindel) + " present(s)"
    elif classmatch not in list_expectedmatrix_header_clade:
        commentary += ""
    if val_poi != '':
        commentary += " ~ " + val_poi
    commentary += '.'
    return commentary


def gen_commentary_rsv_full(classmatch, substitutions, deletions, insertions, nocovsub, nocovins, comment, atypicsub, atypicindel, val_poi):
    if classmatch == 'ININT':
        return 'DSC'
    if nocovsub != []:
        nocovsubpos = [ x.split(':')[0] + ':' + ''.join(filter(str.isdigit, x.split(':')[-1])) for x in nocovsub ]
    else:
        nocovsubpos = []
    commentary = 'Profil des proteines F et G : '
    commentary += ';'.join(substitutions)
    if len(deletions) > 0:
         commentary += '|' + ';'.join(deletions)
    if len(insertions) > 0:
        commentary += '|' + ';'.join(insertions)
    if nocovsub != []:
        if len(list_intersect(nocovsubpos, list_aa_rsvNBS, 'match', 100)) > 0:
            commentary += " ~ Rendu sous reserve compte tenu de(s) position(s) d'interet suivante(s) non couvertes sur les proteines F et G : " + ';'.join(nocovsub)
        else:
            commentary += " ~ Position(s) d'interet suivante(s) non couvertes sur les proteines F et G : " + ';'.join(nocovsub)
    commentary += " ~ " + comment + " (" + expectedmatrix_version[4:] + "." + expectedmatrix_version[2:4] + "." + expectedmatrix_version[0:2] + ")"
    if (atypicsub != [''] or atypicindel != ['']) and classmatch not in list_expectedmatrix_header_clade:
        commentary += " ~ Profil rare :"
        if atypicsub != ['']:
            commentary += " substitution(s) " + ";".join(atypicsub) + " presente(s) sur les proteines F et G"
        if atypicsub != [''] and atypicindel != ['']:
            commentary += " -"
        if atypicindel != ['']:
            commentary += " indel(s) " + ";".join(atypicindel) + " present(s)"
    elif classmatch not in list_expectedmatrix_header_clade:
        commentary += ""
    if val_poi != '':
        commentary += " ~ " + val_poi
    commentary += '.'
    return commentary


def gen_commentary_rsv_simple(classmatch, substitutions, subbaseprofile, deletions, insertions, nocovsub, nocovins, comment, atypicsub, atypicindel, val_poi):
    if classmatch == 'ININT':
        return 'DSC'
    if nocovsub != []:
        nocovsubpos = [ x.split(':')[0] + ':' + ''.join(filter(str.isdigit, x.split(':')[-1])) for x in nocovsub ]
    else:
        nocovsubpos = []
    substitutions_nextclade = [ x.split(':')[0] + ':' + ''.join(filter(str.isdigit, x.split(':')[-1])) for x in substitutions ]
    substitutions = substitutions + [ x for x in subbaseprofile if x.split(':')[0] + ':' + ''.join(filter(str.isdigit, x.split(':')[-1])) not in substitutions_nextclade ]
    substitutions_ordered = [ [ x.split(':')[0], int(re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<3>", x)), re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>\g<4>", x) ] for x in substitutions if x[-1].isalpha() and re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>", x) not in nocovsubpos and re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>", x) in list_aa_rsvNBS ] + [ [ x.split(':')[0], int(re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<3>", x)), re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>\g<2>", x) ] for x in substitutions if x[-1].isdigit() and re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>", x) not in nocovsubpos and re.sub("^([^:]+):([a-zA-Z]+)([0-9]+)([a-zA-Z]*)", "\g<1>:\g<3>", x) in list_aa_rsvNBS ]
    substitutions_ordered = sorted(substitutions_ordered, key=itemgetter(0,1))
    commentary = 'Profil du site de liaison du nirsevimab : '
    commentary += ';'.join([ x[2] for x in substitutions_ordered])
    if len(deletions) > 0:
         commentary += '|' + ';'.join(deletions)
    if len(insertions) > 0:
        commentary += '|' + ';'.join(insertions)
    if nocovsub != []:
        if len(list_intersect(nocovsubpos, list_aa_rsvNBS, 'match', 100)) > 0:
            commentary += " ~ Rendu sous reserve compte tenu de(s) position(s) d'interet suivante(s) non couvertes : " + ';'.join(nocovsub)
        else:
            commentary += " ~ Position(s) d'interet suivante(s) non couvertes : " + ';'.join(nocovsub)
    commentary += " ~ " + comment + " (" + expectedmatrix_version[4:] + "." + expectedmatrix_version[2:4] + "." + expectedmatrix_version[0:2] + ")"
#    if (atypicsub != [''] or atypicindel != ['']) and classmatch not in list_expectedmatrix_header_clade:
#        commentary += " ~ Profil rare :"
#        if atypicsub != ['']:
#            commentary += " substitution(s) " + ";".join(atypicsub) + " presente(s)"
#        if atypicsub != [''] and atypicindel != ['']:
#            commentary += " -"
#        if atypicindel != ['']:
#            commentary += " indel(s) " + ";".join(atypicindel) + " present(s)"
#    elif classmatch not in list_expectedmatrix_header_clade:
#        commentary += ""
#    if val_poi != '':
#        commentary += " ~ " + val_poi
    commentary += '.'
    return commentary

def string_addinfo(string, condition, info):
    if string != condition:
        return string + info
    else:
        return string

parser = argparse.ArgumentParser(description='Generate a validation report from seqmet files')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('--verbose', action='store_true')
debugmode.add_argument('--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.1.5')
parser.add_argument('--runprefix', help='run prefix')
depthdata = parser.add_mutually_exclusive_group()
depthdata.add_argument('--depthfiledir', help='the base', type=path_asfileslist)
depthdata.add_argument('--depthfilelist', help='the base', type=string_asfileslist)
parser.add_argument('--summary', help='the base')
parser.add_argument('--nextclade', help='the base')
parser.add_argument('--pangolin', help='the base')
parser.add_argument('--nextcladeversion', help='the base')
parser.add_argument('--vartable', help='the base')
parser.add_argument('--expectedmatrix', help='the base')
parser.add_argument('--likelymatrix', help='the base')
parser.add_argument('--gff', help='the base')
parser.add_argument('--poitable', help='the base')
parser.add_argument('--mode', help='the base')
parser.add_argument('--varcount_threshold', type=int, default=13, help='the base')
parser.add_argument('--dp_threshold', type=int, default=6, help='the base')
parser.add_argument('--coinf_threshold', type=float, default=2.361, help='the base')
parser.add_argument('--cov_minok', type=int, default=90, help='the base')
parser.add_argument('--cov_maxneg', type=int, default=5, help='the base')
parser.add_argument('--error_plate', default='', help='the base')
parser.add_argument('--error_list', help='the base')
parser.add_argument('--outdir', help='output files to specified dir', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

if args.mode == 'ncov':

    run_date = args.runprefix.split('_')[0]

    datatype = args.expectedmatrix.split('/')[-1].split('_')[1]

    expectedmatrix_version = args.expectedmatrix.split('/')[-1].split('-')[-1].split('.')[0]

    df_var = table_asdf(args.vartable, '\t')
    df_var = df_var.apply(pd.to_numeric, errors='coerce').fillna(df_var)
    df_summary = table_asdf(args.summary, '\t')
    df_nextclade = table_asdf(args.nextclade, '\t')
    df_pangolin = table_asdf(args.pangolin, ',')
    df_expectedmatrix = table_asdf(args.expectedmatrix, '\t')

    list_var = table_aslist(args.vartable, '\t')
    list_var_change = list_var[0][1:]
    list_var_header = [ x[0] for x in list_var ]

    list_poi = table_aslist(args.poitable, '\t')

    list_expectedmatrix = table_aslist(args.expectedmatrix, '\t')
    list_expectedmatrix_change = list_expectedmatrix[0][1:]
    list_expectedmatrix_header = [ x[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_clade = [ x[0].split('_')[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_lineage = [ x[0].split('_')[-1] for x in list_expectedmatrix ]

    list_aarbd = range(333,528)
    list_aarbm = range(438,507)
    list_ntrbd = range(22559,23147)
    list_ntrbm = range(22874,23081)
    list_nts = range(21563,25385)

    list_likelymatrix = table_aslist(args.likelymatrix, '\t')
    list_likelymatrix_change = list_likelymatrix[0][1:]
    list_likelymatrix_header = [ x[0] for x in list_likelymatrix ]
    list_likelymatrix_header_clade = [ x[0].split('_')[0] for x in list_likelymatrix ]
    list_likelymatrix_header_lineage = [ x[0].split('_')[-1] for x in list_likelymatrix ]

    list_profile_base = list_expectedmatrix[0][1:]
    
    list_likely_base = list_likelymatrix[0][1:]
    
    list_gffposlist = gff_asposlist(args.gff, 'gene', ['CDS'])

    list_errorplate = args.error_plate.split(',')

    list_errorlist = open(args.error_list, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n') if args.error_list else []

    df_summary.rename(columns=lambda x: 'summary_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_nextclade.rename(columns=lambda x: 'nextclade_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_pangolin.rename(columns=lambda x: 'pangolin_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_validation = pd.concat([df_summary, df_nextclade, df_pangolin], axis=1)

    df_validation['val_dp'] = df_validation.apply(lambda x : seek_dp(x.summary_vcf_dpcount, '/', args.dp_threshold) , axis=1)

    df_validation['val_varcount'] = df_validation[['summary_vcf_af10-20count', 'summary_vcf_af20-50count']].apply(lambda x: pd.to_numeric(x, errors='coerce')).sum(axis = 1, skipna = False, min_count=2).apply(lambda x : seek_varcount(x, '//', args.varcount_threshold))

    df_validation['stat_tableindex'] = df_validation.apply(lambda x : list_var_header.index(x.name) if x.name in list_var_header else 0, axis=1)

    df_validation['val_poi'] = df_validation.apply(lambda x : seek_poi(list_var_change, list_var[x.stat_tableindex][1:], list_poi) , axis=1)

    df_validation['val_plate'] = df_validation.apply(lambda x : x.name.split('-')[0] if '-' in x.name else 'Pl000', axis=1)

    df_validation['val_glims'] = df_validation.apply(lambda x : x.name.split('-')[1].split('_')[0] if '-' in x.name else x.name.split('_')[0], axis=1)

    df_validation['val_platewell'] = df_validation.apply(lambda x : x.name.split('_')[-1], axis=1)

    df_validation['stat_missingabspos'] = df_validation.apply(lambda x : list_position(str(x.nextclade_missing), '-', ',') , axis=1)

#    df_validation['stat_missingntpos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'nt') , axis=1)

    df_validation['stat_missingaapos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'aa') , axis=1)

    df_validation['val_rbmcoverage'] = df_validation.apply(lambda x : list_intersect(x.stat_missingabspos, list_ntrbm, 'percent', 100) , axis=1)

    df_validation['val_spikecoverage'] = df_validation.apply(lambda x : list_intersect(x.stat_missingabspos, list_nts, 'percent', 100) , axis=1)

    df_validation['val_clade'] = df_validation.apply(lambda x : validate_column_ncov(x.nextclade_clade_legacy, x.summary_consensus_perccoverage, x.val_spikecoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['val_lineage'] = df_validation.apply(lambda x : validate_column_ncov(x.pangolin_lineage, x.summary_consensus_perccoverage, x.val_spikecoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['val_classmatch'] = df_validation.apply(lambda x : match_expectedmatrix_ncov(x.val_lineage, x.val_clade, list_expectedmatrix_header, "match") , axis=1)

    df_validation['val_classindex'] = df_validation.apply(lambda x : match_expectedmatrix_ncov(x.val_lineage, x.val_clade, list_expectedmatrix_header, "index") , axis=1)

    df_validation['val_classcomment'] = df_validation.apply(lambda x : match_expectedmatrix_ncov(x.val_lineage, x.val_clade, list_expectedmatrix_header, "comment") , axis=1)

    df_validation['val_insertions'] = df_validation.apply(lambda x : seek_insertion_ncov(str(x.nextclade_insertions)) , axis=1)

    df_validation['val_expectedprofile'] = df_validation.apply(lambda x : seek_profile(list_expectedmatrix_change, list_expectedmatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_likelyprofile'] = df_validation.apply(lambda x : seek_profile(list_likelymatrix_change, list_likelymatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_nocovsub'] = df_validation.apply(lambda x : seek_nocovsub(list_expectedmatrix_change, [ y for y in x.stat_missingaapos if y[0:2] == 'S:']) , axis=1)

    df_validation['val_nocovins'] = df_validation.apply(lambda x : seek_nocovins(list_expectedmatrix_change, x.stat_missingabspos) , axis=1)

    df_validation['val_missingsub'] = df_validation.apply(lambda x : list_intersect(x.val_expectedprofile, [y for y in (str(x.nextclade_aaSubstitutions) + "," + str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',')], 'missing', 100) , axis=1)

    df_validation['val_atypicsub'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaSubstitutions)).split(',') if y[0:2] == 'S:' and (int(''.join(filter(str.isdigit, y)))) in list_aarbd and y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_atypicindel'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if (y[0:2] == 'S:' or y[0:3] == 'Ins') and y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_avisbio'] = df_validation.apply(lambda x : seek_issue_ncov(x.val_glims, x.val_clade, x.summary_qc_seqcontrol, x.val_spikecoverage, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_nocovsub, x.val_nocovins, x.val_missingsub, x.nextclade_qc_privateMutations_status, x.summary_vcf_coinf01score, x.summary_vcf_coinf02match, x.summary_vcf_coinf02count, x.summary_vcf_coinf02score, x.val_poi, 'x.summary_bam_verif', args.coinf_threshold) , axis=1)

    df_validation['val_result'] = df_validation.apply(lambda x : validate_result(x.val_clade, x.val_dp, x.val_varcount, x.summary_qc_seqcontrol, x.val_avisbio) , axis=1)

    df_validation['val_commentary'] = df_validation.apply(lambda x : gen_commentary_ncov(x.val_classmatch, [y for y in str(x.nextclade_aaSubstitutions).split(',') if y not in ['', 'nan'] and y.split(':')[0] == 'S'], [y for y in str(x.nextclade_aaDeletions).split(',') if y not in ['', 'nan'] and y.split(':')[0] == 'S'], [y for y in str(x.val_insertions).split(',') if y not in ['', 'nan'] and int(y.split(':')[-1].rstrip('ACGTN')) in list_nts], x.val_nocovsub, x.val_nocovins, x.val_classcomment, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_poi), axis=1)

    df_validation['val_error'] = df_validation.apply(lambda x : 'ERROR' if x.val_plate in list_errorplate or x.name in list_errorlist else "OK", axis=1)

    df_validation['stat_simpleprofile'] = df_validation.apply(lambda x : x.val_clade + ' - ' + x.val_atypicsub , axis=1)

    df_validation['stat_pangolin'] = df_validation.apply(lambda x : string_addinfo(x.val_lineage, 'ININT', ' (' + str(x.pangolin_version) + ')') , axis=1)

    df_var_major = filter_dfrownames(df_var, ['|MAJOR'], 'keep')
    df_var_major_indel = filter_dfrownames(df_var_major, ['|DEL|', '|INS|', '|COMPLEX|'], 'keep')

    df_validation_sample = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'drop')
    df_validation_control = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_nt = filter_dfrownames(df_validation, ['NT', 'PCR', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_20a = filter_dfrownames(df_validation, ['TposPl'], 'keep')
    df_validation_20h = filter_dfrownames(df_validation, ['TposV2'], 'keep')
    df_validation_20j = filter_dfrownames(df_validation, ['TposV3'], 'keep')
    df_validation_19b = filter_dfrownames(df_validation, ['Tpos19B'], 'keep')

    df_var_major_indel[df_var_major_indel > 0].count(axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_indel.tsv', sep='\t', index_label='variant', header=["count"])

    df_validation_control[['summary_consensus_perccoverage', 'nextclade_clade_legacy', 'nextclade_totalAminoacidSubstitutions']].replace(r'^\s*$', np.nan, regex=True).fillna('ININT').to_csv(args.outdir + '/' + run_date + '_' + datatype + '_control.tsv', sep='\t', index_label='control')

    df_validation[['val_result', 'val_lineage', 'val_commentary', 'val_avisbio', 'val_nocovsub', 'val_nocovins', 'val_missingsub', 'val_atypicsub', 'val_atypicindel', 'val_spikecoverage', 'summary_bam_meandepth', 'summary_consensus_perccoverage', 'summary_vcf_dpcount', 'summary_vcf_af05-50count', 'summary_vcf_coinf01match', 'summary_vcf_coinf01count', 'summary_vcf_coinf01score', 'summary_vcf_coinf02match', 'summary_vcf_coinf02count', 'summary_vcf_coinf02score', 'summary_run_id']].replace(r'^\s*$', np.nan, regex=True).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_shortval.tsv', sep='\t', index_label='sample_id')

    df_validation_sample['stat_simpleprofile'].value_counts().sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_runsummary.tsv', sep='\t', index_label='Clade')

    df_metrics = pd.concat([
        df_validation_sample[df_validation_sample['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') < args.cov_minok]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_dp'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_varcount'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample[df_validation_sample['summary_qc_seqcontrol'].isin(['FAILED'])]['val_plate'].value_counts(),
        df_validation_nt[df_validation_nt['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg]['val_plate'].value_counts(),
        df_validation_20a.loc[(df_validation_20a['nextclade_clade_legacy'].isin(['20A'])) & (df_validation_20a['nextclade_totalAminoacidSubstitutions'].apply(pd.to_numeric, errors='coerce') <= 6.0)]['val_plate'].value_counts().add(df_validation_20h[df_validation_20h['nextclade_clade_legacy'].isin(['20H (Beta, V2)'])]['val_plate'].value_counts(), fill_value=0).add(df_validation_20j[df_validation_20j['nextclade_clade_legacy'].isin(['20J (Gamma, V3)'])]['val_plate'].value_counts(), fill_value=0).add(df_validation_19b[df_validation_19b['nextclade_clade_legacy'].isin(['19B'])]['val_plate'].value_counts(), fill_value=0)],axis=1).fillna(0).astype(int)
    df_metrics.columns = ['sample_cov<=' + str(args.cov_minok), 'sample_dp>=' + str(args.dp_threshold), 'sample_varcount>=' + str(args.varcount_threshold), 'sample_seq_failed', 'tneg_cov>=' + str(args.cov_maxneg), 'tpos_ok']
    df_metrics.sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_metric.tsv', sep='\t', index_label='plate')

    pd.concat([df_validation[df_validation.filter(regex='summary_').columns], df_validation[df_validation.filter(regex='val_').columns], df_validation[df_validation.filter(regex='nextclade_').columns], df_validation[df_validation.filter(regex='pangolin_').columns]], axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_validation.tsv', sep='\t', index_label='sample_id')

    df_export = pd.DataFrame(df_validation, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "SeqArt"
    df_export['Target_1_cq'] = df_validation.apply(lambda x : x.val_result.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export['Target_2_cq'] = df_validation.apply(lambda x : x.stat_pangolin.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.stat_pangolin.replace(',',';'), axis=1)
    df_export['Target_3_cq'] = df_validation.apply(lambda x : x.val_commentary.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_commentary.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False)



if args.mode == 'fluabv':

    run_date = args.runprefix.split('_')[0]

    datatype = args.expectedmatrix.split('/')[-1].split('_')[1]

    expectedmatrix_version = args.expectedmatrix.split('/')[-1].split('-')[-1].split('.')[0]

    df_var = table_asdf(args.vartable, '\t')
    df_var = df_var.apply(pd.to_numeric, errors='coerce').fillna(df_var)
    df_summary = table_asdf(args.summary, '\t')
    df_nextclade = table_asdf(args.nextclade, '\t')
    df_expectedmatrix = table_asdf(args.expectedmatrix, '\t')

    list_var = table_aslist(args.vartable, '\t')
    list_var_change = list_var[0][1:]
    list_var_header = [ x[0] for x in list_var ]

    list_poi = table_aslist(args.poitable, '\t')

    list_expectedmatrix = table_aslist(args.expectedmatrix, '\t')
    list_expectedmatrix_change = list_expectedmatrix[0][1:]
    list_expectedmatrix_header = [ x[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_clade = [ x[0].split('_')[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_expectedmatrix ]

    list_likelymatrix = table_aslist(args.likelymatrix, '\t')
    list_likelymatrix_change = list_likelymatrix[0][1:]
    list_likelymatrix_header = [ x[0] for x in list_likelymatrix ]
    list_likelymatrix_header_clade = [ x[0].split('_')[0] for x in list_likelymatrix ]
    list_likelymatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_likelymatrix ]

    list_profile_base = list_expectedmatrix[0][1:]
    
    list_likely_base = list_likelymatrix[0][1:]
    
    list_gffposlist = gff_asposlist(args.gff, 'gene', ['CDS'])

    list_errorplate = args.error_plate.split(',')

    list_errorlist = open(args.error_list, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n') if args.error_list else []

    df_summary.rename(columns=lambda x: 'summary_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_nextclade.rename(columns=lambda x: 'nextclade_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_validation = pd.concat([df_summary, df_nextclade], axis=1)

    df_validation['val_dp'] = df_validation.apply(lambda x : seek_dp(x.summary_vcf_dpcount, '/', args.dp_threshold) , axis=1)

    df_validation['val_varcount'] = df_validation[['summary_vcf_af10-20count', 'summary_vcf_af20-50count']].apply(lambda x: pd.to_numeric(x, errors='coerce')).sum(axis = 1, skipna = False, min_count=2).apply(lambda x : seek_varcount(x, '//', args.varcount_threshold))

    df_validation['stat_tableindex'] = df_validation.apply(lambda x : list_var_header.index(x.name) if x.name in list_var_header else 0, axis=1)

    df_validation['val_poi'] = df_validation.apply(lambda x : seek_poi(list_var_change, list_var[x.stat_tableindex][1:], list_poi) , axis=1)

    df_validation['val_plate'] = df_validation.apply(lambda x : x.name.split('-')[0] if '-' in x.name else 'Pl000', axis=1)

    df_validation['val_glims'] = df_validation.apply(lambda x : x.name.split('-')[1].split('_')[0] if '-' in x.name else x.name.split('_')[0], axis=1)

    df_validation['val_platewell'] = df_validation.apply(lambda x : x.name.split('_')[-1] , axis=1)

    df_validation['val_clade'] = df_validation.apply(lambda x : validate_column_fluabv(x.nextclade_clade, x.summary_consensus_perccoverage_S4, x.summary_consensus_perccoverage_S6, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['stat_missingabspos'] = df_validation.apply(lambda x : list_position_fluabv(str(x.nextclade_missing), '-', ',', str(x.nextclade_insertions)) , axis=1)

#    df_validation['stat_missingntpos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'nt') , axis=1)

    df_validation['stat_missingaapos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'aa') , axis=1)

    df_validation['val_classmatch'] = df_validation.apply(lambda x : match_expectedmatrix_fluabv(x.val_clade, '', list_expectedmatrix_header, "match") , axis=1)

    df_validation['val_classindex'] = df_validation.apply(lambda x : match_expectedmatrix_fluabv(x.val_clade, '', list_expectedmatrix_header, "index") , axis=1)

    df_validation['val_classcomment'] = df_validation.apply(lambda x : match_expectedmatrix_fluabv(x.val_clade, '', list_expectedmatrix_header, "comment") , axis=1)

    df_validation['val_insertions'] = df_validation.apply(lambda x : seek_insertion(str(x.nextclade_insertions)) , axis=1)

    df_validation['val_expectedprofile'] = df_validation.apply(lambda x : seek_profile(list_expectedmatrix_change, list_expectedmatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_likelyprofile'] = df_validation.apply(lambda x : seek_profile(list_likelymatrix_change, list_likelymatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_nocovsub'] = df_validation.apply(lambda x : seek_nocovsub(list_expectedmatrix_change, x.stat_missingaapos) , axis=1)

    df_validation['val_nocovins'] = df_validation.apply(lambda x : seek_nocovins(list_expectedmatrix_change, x.stat_missingabspos) , axis=1)

    df_validation['val_missingsub'] = df_validation.apply(lambda x : list_intersect(x.val_expectedprofile, [y for y in (str(x.nextclade_aaSubstitutions) + "," + str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if y not in ['', 'nan']], 'missing', 100) , axis=1)

    df_validation['val_atypicsub'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaSubstitutions)).split(',') if y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_atypicindel'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_avisbio'] = df_validation.apply(lambda x : seek_issue_fluabv(x.val_glims, x.val_clade, x.summary_qc_seqcontrol, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_nocovsub, x.val_nocovins, x.val_missingsub, x.nextclade_qc_privateMutations_status, x.summary_vcf_coinf01score, x.summary_vcf_coinf02match, x.summary_vcf_coinf02count, x.summary_vcf_coinf02score, x.val_poi, x.summary_bam_verif, args.coinf_threshold) , axis=1)

    df_validation['val_result'] = df_validation.apply(lambda x : validate_result_fluabv(x.val_clade, x.val_dp, x.val_varcount, x.summary_qc_seqcontrol, x.val_avisbio) , axis=1)

    df_validation['val_commentary'] = df_validation.apply(lambda x : gen_commentary_fluabv(x.val_classmatch, [y for y in str(x.nextclade_aaSubstitutions).split(',') if y not in ['', 'nan'] and y[:-1] in list_profile_base], [y for y in str(x.nextclade_aaDeletions).split(',') if y not in ['', 'nan'] and y[:-1] in list_profile_base], [y for y in str(x.val_insertions).split(',') if y not in ['', 'nan'] and y.rstrip('ACGTN') in list_profile_base], x.val_nocovsub, x.val_nocovins, x.val_classcomment, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_poi), axis=1)

    df_validation['val_error'] = df_validation.apply(lambda x : 'ERROR' if x.val_plate in list_errorplate or x.name in list_errorlist else "OK", axis=1)

    df_validation['stat_simpleprofile'] = df_validation.apply(lambda x : x.val_clade + ' - ' + x.val_atypicsub , axis=1)

    df_var_major = filter_dfrownames(df_var, ['|MAJOR'], 'keep')
    df_var_major_indel = filter_dfrownames(df_var_major, ['|DEL|', '|INS|', '|COMPLEX|'], 'keep')

    df_validation_sample = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'drop')
    df_validation_control = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_nt = filter_dfrownames(df_validation, ['NT', 'PCR', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')

    df_var_major_indel[df_var_major_indel > 0].count(axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_indel.tsv', sep='\t', index_label='variant', header=["count"])

    df_validation_control[['summary_consensus_perccoverage_S4', 'summary_consensus_perccoverage_S6', 'nextclade_clade', 'nextclade_totalAminoacidSubstitutions']].replace(r'^\s*$', np.nan, regex=True).fillna('ININT').to_csv(args.outdir + '/' + run_date + '_' + datatype + '_control.tsv', sep='\t', index_label='control')

    df_validation[['val_result', 'val_commentary', 'val_avisbio', 'val_nocovsub', 'val_nocovins', 'val_missingsub', 'val_atypicsub', 'val_atypicindel', 'summary_bam_meandepth_S1', 'summary_bam_meandepth_S2', 'summary_bam_meandepth_S3', 'summary_bam_meandepth_S4', 'summary_bam_meandepth_S5', 'summary_bam_meandepth_S6', 'summary_bam_meandepth_S7', 'summary_bam_meandepth_S8', 'summary_consensus_perccoverage_S1', 'summary_consensus_perccoverage_S2', 'summary_consensus_perccoverage_S3', 'summary_consensus_perccoverage_S4', 'summary_consensus_perccoverage_S5', 'summary_consensus_perccoverage_S6', 'summary_consensus_perccoverage_S7', 'summary_consensus_perccoverage_S8', 'summary_vcf_dpcount', 'summary_vcf_af05-50count', 'summary_vcf_coinf01match', 'summary_vcf_coinf01count', 'summary_vcf_coinf01score', 'summary_vcf_coinf02match', 'summary_vcf_coinf02count', 'summary_vcf_coinf02score', 'summary_run_id']].replace(r'^\s*$', np.nan, regex=True).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_shortval.tsv', sep='\t', index_label='sample_id')

    df_validation_sample['stat_simpleprofile'].value_counts().sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_runsummary.tsv', sep='\t', index_label='Clade')

    df_metrics = pd.concat([
        df_validation_sample.loc[(df_validation_sample['summary_consensus_perccoverage_S4'].apply(pd.to_numeric, errors='coerce') < args.cov_minok) | (df_validation_sample['summary_consensus_perccoverage_S6'].apply(pd.to_numeric, errors='coerce') < args.cov_minok)]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_dp'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_varcount'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample[df_validation_sample['summary_qc_seqcontrol'].isin(['FAILED'])]['val_plate'].value_counts(),
        df_validation_nt.loc[(df_validation_nt['summary_consensus_perccoverage_S4'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg) | (df_validation_nt['summary_consensus_perccoverage_S6'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg)]['val_plate'].value_counts()],axis=1).fillna(0).astype(int)
    df_metrics.columns = ['sample_cov<=' + str(args.cov_minok), 'sample_dp>=' + str(args.dp_threshold), 'sample_varcount>=' + str(args.varcount_threshold), 'sample_seq_failed', 'tneg_cov>=' + str(args.cov_maxneg)]
    df_metrics.sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_metric.tsv', sep='\t', index_label='plate')

    pd.concat([df_validation[df_validation.filter(regex='summary_').columns], df_validation[df_validation.filter(regex='val_').columns], df_validation[df_validation.filter(regex='nextclade_').columns]], axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_validation.tsv', sep='\t', index_label='sample_id')
    #df_validation_sample_notna = df_validation_sample[~df_validation_sample['summary_consensus_perccoverage_S4'].isin(['NA'])]

    df_export = pd.DataFrame(df_validation_sample, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "GABILL"
    df_export['AssayResult'] = df_validation_sample.apply(lambda x : x.val_result.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False)

    df_export = pd.DataFrame(df_validation_sample, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "COMGRAB"
    df_export['AssayResult'] = df_validation_sample.apply(lambda x : x.val_commentary.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False, mode='a', header=False)



if args.mode == 'hmpxv':

    run_date = args.runprefix.split('_')[0]

    datatype = args.expectedmatrix.split('/')[-1].split('_')[1]

    expectedmatrix_version = args.expectedmatrix.split('/')[-1].split('-')[-1].split('.')[0]

    df_var = table_asdf(args.vartable, '\t')
    df_var = df_var.apply(pd.to_numeric, errors='coerce').fillna(df_var)
    df_summary = table_asdf(args.summary, '\t')
    df_nextclade = table_asdf(args.nextclade, '\t')
    df_expectedmatrix = table_asdf(args.expectedmatrix, '\t')

    list_var = table_aslist(args.vartable, '\t')
    list_var_change = list_var[0][1:]
    list_var_header = [ x[0] for x in list_var ]

    list_poi = table_aslist(args.poitable, '\t')

    list_expectedmatrix = table_aslist(args.expectedmatrix, '\t')
    list_expectedmatrix_change = list_expectedmatrix[0][1:]
    list_expectedmatrix_header = [ x[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_clade = [ x[0].split('_')[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_expectedmatrix ]

    list_likelymatrix = table_aslist(args.likelymatrix, '\t')
    list_likelymatrix_change = list_likelymatrix[0][1:]
    list_likelymatrix_header = [ x[0] for x in list_likelymatrix ]
    list_likelymatrix_header_clade = [ x[0].split('_')[0] for x in list_likelymatrix ]
    list_likelymatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_likelymatrix ]

    list_profile_base = list_expectedmatrix[0][1:]
    
    list_likely_base = list_likelymatrix[0][1:]
    
    list_gffposlist = gff_asposlist(args.gff, 'gene', ['CDS'])

    list_errorplate = args.error_plate.split(',')

    list_errorlist = open(args.error_list, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n') if args.error_list else []

    df_summary.rename(columns=lambda x: 'summary_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_nextclade.rename(columns=lambda x: 'nextclade_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_validation = pd.concat([df_summary, df_nextclade], axis=1)

    df_validation['val_dp'] = df_validation.apply(lambda x : seek_dp(x.summary_vcf_dpcount, '/', args.dp_threshold) , axis=1)

    df_validation['val_varcount'] = df_validation[['summary_vcf_af10-20count', 'summary_vcf_af20-50count']].apply(lambda x: pd.to_numeric(x, errors='coerce')).sum(axis = 1, skipna = False, min_count=2).apply(lambda x : seek_varcount(x, '//', args.varcount_threshold))

    df_validation['stat_tableindex'] = df_validation.apply(lambda x : list_var_header.index(x.name) if x.name in list_var_header else 0, axis=1)

    df_validation['val_poi'] = df_validation.apply(lambda x : seek_poi(list_var_change, list_var[x.stat_tableindex][1:], list_poi) , axis=1)

    df_validation['val_plate'] = df_validation.apply(lambda x : x.name.split('-')[0] if '-' in x.name else 'Pl000', axis=1)

    df_validation['val_glims'] = df_validation.apply(lambda x : x.name.split('-')[1].split('_')[0] if '-' in x.name else x.name.split('_')[0], axis=1)

    df_validation['val_platewell'] = df_validation.apply(lambda x : x.name.split('_')[-1] , axis=1)

    df_validation['val_clade'] = df_validation.apply(lambda x : validate_column_hmpxv(x.nextclade_clade, x.summary_consensus_perccoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['val_Gclade'] = df_validation.apply(lambda x : validate_column_hmpxv(x.nextclade_G_clade, x.summary_consensus_perccoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['stat_missingabspos'] = df_validation.apply(lambda x : list_position(str(x.nextclade_missing), '-', ',') , axis=1)

#    df_validation['stat_missingntpos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'nt') , axis=1)

    df_validation['stat_missingaapos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'aa') , axis=1)

    df_validation['val_classmatch'] = df_validation.apply(lambda x : match_expectedmatrix_hmpxv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "match") , axis=1)

    df_validation['val_classindex'] = df_validation.apply(lambda x : match_expectedmatrix_hmpxv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "index") , axis=1)

    df_validation['val_classcomment'] = df_validation.apply(lambda x : match_expectedmatrix_hmpxv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "comment") , axis=1)

    df_validation['val_insertions'] = df_validation.apply(lambda x : seek_insertion(str(x.nextclade_insertions)) , axis=1)

    df_validation['val_expectedprofile'] = df_validation.apply(lambda x : seek_profile(list_expectedmatrix_change, list_expectedmatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_likelyprofile'] = df_validation.apply(lambda x : seek_profile(list_likelymatrix_change, list_likelymatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_nocovsub'] = df_validation.apply(lambda x : seek_nocovsub(list_expectedmatrix_change, x.stat_missingaapos) , axis=1)

    df_validation['val_nocovins'] = df_validation.apply(lambda x : seek_nocovins(list_expectedmatrix_change, x.stat_missingabspos) , axis=1)

    df_validation['val_missingsub'] = df_validation.apply(lambda x : list_intersect(x.val_expectedprofile, [y for y in (str(x.nextclade_aaSubstitutions) + "," + str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if y not in ['', 'nan']], 'missing', 100) , axis=1)

    df_validation['val_atypicsub'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaSubstitutions)).split(',') if y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_atypicindel'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_avisbio'] = df_validation.apply(lambda x : seek_issue_hmpxv(x.val_glims, x.val_clade, x.summary_qc_seqcontrol, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_nocovsub, x.val_nocovins, x.val_missingsub, x.nextclade_qc_privateMutations_status, x.summary_vcf_coinf01score, x.summary_vcf_coinf02match, x.summary_vcf_coinf02count, x.summary_vcf_coinf02score, x.val_poi, 'x.summary_bam_verif', args.coinf_threshold) , axis=1)

    df_validation['val_result'] = df_validation.apply(lambda x : validate_result_hmpxv(x.val_clade, x.val_Gclade, x.val_dp, x.val_varcount, x.summary_qc_seqcontrol, x.val_avisbio) , axis=1)

    df_validation['val_commentary'] = df_validation.apply(lambda x : gen_commentary_hmpxv(x.val_classmatch, [y for y in str(x.nextclade_aaSubstitutions).split(',') if y not in ['', 'nan'] and y[:-1] in list_profile_base], [y for y in str(x.nextclade_aaDeletions).split(',') if y not in ['', 'nan'] and y[:-1] in list_profile_base], [y for y in str(x.val_insertions).split(',') if y not in ['', 'nan'] and y.rstrip('ACGTN') in list_profile_base], x.val_nocovsub, x.val_nocovins, x.val_classcomment, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_poi), axis=1)

    df_validation['val_error'] = df_validation.apply(lambda x : 'ERROR' if x.val_plate in list_errorplate or x.name in list_errorlist else "OK", axis=1)

    df_validation['stat_simpleprofile'] = df_validation.apply(lambda x : x.val_clade + ' - ' + x.val_atypicsub , axis=1)

    df_var_major = filter_dfrownames(df_var, ['|MAJOR'], 'keep')
    df_var_major_indel = filter_dfrownames(df_var_major, ['|DEL|', '|INS|', '|COMPLEX|'], 'keep')

    df_validation_sample = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'drop')
    df_validation_control = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_nt = filter_dfrownames(df_validation, ['NT', 'PCR', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_tpos = filter_dfrownames(df_validation, ['Tpos'], 'keep')

    df_var_major_indel[df_var_major_indel > 0].count(axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_indel.tsv', sep='\t', index_label='variant', header=["count"])

    df_validation_control[['summary_consensus_perccoverage', 'nextclade_clade', 'nextclade_totalAminoacidSubstitutions']].replace(r'^\s*$', np.nan, regex=True).fillna('ININT').to_csv(args.outdir + '/' + run_date + '_' + datatype + '_control.tsv', sep='\t', index_label='control')

    df_validation[['val_result', 'val_Gclade', 'val_commentary', 'val_avisbio', 'val_nocovsub', 'val_nocovins', 'val_missingsub', 'val_atypicsub', 'val_atypicindel', 'summary_bam_meandepth', 'summary_consensus_perccoverage', 'summary_vcf_dpcount', 'summary_vcf_af05-50count', 'summary_vcf_coinf01match', 'summary_vcf_coinf01count', 'summary_vcf_coinf01score', 'summary_vcf_coinf02match', 'summary_vcf_coinf02count', 'summary_vcf_coinf02score', 'summary_run_id']].replace(r'^\s*$', np.nan, regex=True).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_shortval.tsv', sep='\t', index_label='sample_id')

    df_validation_sample['stat_simpleprofile'].value_counts().sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_runsummary.tsv', sep='\t', index_label='Clade')

    df_metrics = pd.concat([
        df_validation_sample[df_validation_sample['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') < args.cov_minok]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_dp'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_varcount'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample[df_validation_sample['summary_qc_seqcontrol'].isin(['FAILED'])]['val_plate'].value_counts(),
        df_validation_nt.loc[(df_validation_nt['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg)]['val_plate'].value_counts()],axis=1).fillna(0).astype(int)
    df_metrics.columns = ['sample_cov<=' + str(args.cov_minok), 'sample_dp>=' + str(args.dp_threshold), 'sample_varcount>=' + str(args.varcount_threshold), 'sample_seq_failed', 'tneg_cov>=' + str(args.cov_maxneg)]
    df_metrics.sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_metric.tsv', sep='\t', index_label='plate')

    pd.concat([df_validation[df_validation.filter(regex='summary_').columns], df_validation[df_validation.filter(regex='val_').columns], df_validation[df_validation.filter(regex='nextclade_').columns]], axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_validation.tsv', sep='\t', index_label='sample_id')

    df_export = pd.DataFrame(df_validation, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "MPXVMETA"
    df_export['AssayResult'] = df_validation.apply(lambda x : x.val_result.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False)


if args.mode == 'rsv':

    run_date = args.runprefix.split('_')[0]

    datatype = args.expectedmatrix.split('/')[-1].split('_')[1]

    expectedmatrix_version = args.expectedmatrix.split('/')[-1].split('-')[-1].split('.')[0]

    df_var = table_asdf(args.vartable, '\t')
    df_var = df_var.apply(pd.to_numeric, errors='coerce').fillna(df_var)
    df_summary = table_asdf(args.summary, '\t')
    df_nextclade = table_asdf(args.nextclade, '\t')
    df_expectedmatrix = table_asdf(args.expectedmatrix, '\t')

    list_var = table_aslist(args.vartable, '\t')
    list_var_change = list_var[0][1:]
    list_var_header = [ x[0] for x in list_var ]

    list_poi = table_aslist(args.poitable, '\t')

    list_expectedmatrix = table_aslist(args.expectedmatrix, '\t')
    list_expectedmatrix_change = list_expectedmatrix[0][1:]
    list_expectedmatrix_header = [ x[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_clade = [ x[0].split('_')[0] for x in list_expectedmatrix ]
    list_expectedmatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_expectedmatrix ]

    list_likelymatrix = table_aslist(args.likelymatrix, '\t')
    list_likelymatrix_change = list_likelymatrix[0][1:]
    list_likelymatrix_header = [ x[0] for x in list_likelymatrix ]
    list_likelymatrix_header_clade = [ x[0].split('_')[0] for x in list_likelymatrix ]
    list_likelymatrix_header_Gclade = [ x[0].split('_')[-1] for x in list_likelymatrix ]

    list_profile_base = list_expectedmatrix[0][1:]
    
    list_likely_base = list_likelymatrix[0][1:]
    
    list_gffposlist = gff_asposlist(args.gff, 'gene', ['CDS'])
    list_gffposlist_NS1 = list_gffposlist[[ x[0] for x in list_gffposlist ].index('NS1')]
    list_gffposlist_NS2 = list_gffposlist[[ x[0] for x in list_gffposlist ].index('NS2')]
    list_gffposlist_N = list_gffposlist[[ x[0] for x in list_gffposlist ].index('N')]
    list_gffposlist_P = list_gffposlist[[ x[0] for x in list_gffposlist ].index('P')]
    list_gffposlist_M = list_gffposlist[[ x[0] for x in list_gffposlist ].index('M')]
    list_gffposlist_SH = list_gffposlist[[ x[0] for x in list_gffposlist ].index('SH')]
    list_gffposlist_G = list_gffposlist[[ x[0] for x in list_gffposlist ].index('G')]
    list_gffposlist_F = list_gffposlist[[ x[0] for x in list_gffposlist ].index('F')]
    list_gffposlist_M21 = list_gffposlist[[ x[0] for x in list_gffposlist ].index('M2-1')]
    list_gffposlist_M22 = list_gffposlist[[ x[0] for x in list_gffposlist ].index('M2-2')]
    list_gffposlist_L = list_gffposlist[[ x[0] for x in list_gffposlist ].index('L')]

    list_aa_rsvG = list_gffposlist_G[2]
    list_nt_rsvG = list_gffposlist_G[1]
    list_aa_rsvF = list_gffposlist_F[2]
    list_nt_rsvF = list_gffposlist_F[1]
    list_aa_rsvNBS = ['F:62', 'F:63', 'F:64', 'F:65', 'F:66', 'F:67', 'F:68', 'F:69', 'F:196', 'F:197', 'F:198', 'F:199', 'F:200', 'F:201', 'F:202', 'F:203', 'F:204', 'F:205', 'F:206', 'F:207', 'F:208', 'F:209', 'F:210', 'F:211', 'F:212']
    list_nt_cds = list_gffposlist_NS1[1]+list_gffposlist_NS2[1]+list_gffposlist_N[1]+list_gffposlist_P[1]+list_gffposlist_M[1]+list_gffposlist_SH[1]+list_gffposlist_G[1]+list_gffposlist_F[1]+list_gffposlist_M21[1]+list_gffposlist_M22[1]+list_gffposlist_L[1]
    list_target_gene = ['G', 'F']
#    list_target_gene = [ x[0] for x in list_gffposlist ]

    list_errorplate = args.error_plate.split(',')

    list_errorlist = open(args.error_list, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n') if args.error_list else []

    df_summary.rename(columns=lambda x: 'summary_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_nextclade.rename(columns=lambda x: 'nextclade_' + x.replace('%','perc').replace('.','_'), inplace=True)
    df_validation = pd.concat([df_summary, df_nextclade], axis=1)

    df_validation['val_dp'] = df_validation.apply(lambda x : seek_dp(x.summary_vcf_dpcount, '/', args.dp_threshold) , axis=1)

    df_validation['val_varcount'] = df_validation[['summary_vcf_af10-20count', 'summary_vcf_af20-50count']].apply(lambda x: pd.to_numeric(x, errors='coerce')).sum(axis = 1, skipna = False, min_count=2).apply(lambda x : seek_varcount(x, '//', args.varcount_threshold))

    df_validation['stat_tableindex'] = df_validation.apply(lambda x : list_var_header.index(x.name) if x.name in list_var_header else 0, axis=1)

    df_validation['val_poi'] = df_validation.apply(lambda x : seek_poi(list_var_change, list_var[x.stat_tableindex][1:], list_poi) , axis=1)

    df_validation['val_plate'] = df_validation.apply(lambda x : x.name.split('-')[0] if '-' in x.name else 'Pl000', axis=1)

    df_validation['val_glims'] = df_validation.apply(lambda x : x.name.split('-')[1].split('_')[0] if '-' in x.name else x.name.split('_')[0], axis=1)

    df_validation['val_platewell'] = df_validation.apply(lambda x : x.name.split('_')[-1] , axis=1)

    df_validation['stat_missingabspos'] = df_validation.apply(lambda x : list_position(str(x.nextclade_missing), '-', ',') , axis=1)

    #df_validation['stat_missingntpos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'nt') , axis=1)

    df_validation['stat_missingaapos'] = df_validation.apply(lambda x : seek_missingpos(x.stat_missingabspos, list_gffposlist, 'aa') , axis=1)

    df_validation['val_Fcoverage'] = df_validation.apply(lambda x : list_intersect(x.stat_missingabspos, list_nt_rsvF, 'percent', 100) , axis=1)

    df_validation['val_Gcoverage'] = df_validation.apply(lambda x : list_intersect(x.stat_missingabspos, list_nt_rsvG, 'percent', 100) , axis=1)

    df_validation['val_clade'] = df_validation.apply(lambda x : validate_column_rsv(x.nextclade_clade, x.summary_consensus_perccoverage, x.val_Fcoverage, x.val_Gcoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['val_Gclade'] = df_validation.apply(lambda x : validate_column_rsv(x.nextclade_G_clade, x.summary_consensus_perccoverage, x.val_Fcoverage, x.val_Gcoverage, ';', args.cov_minok, x.nextclade_qc_overallStatus) , axis=1)

    df_validation['val_classmatch'] = df_validation.apply(lambda x : match_expectedmatrix_rsv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "match") , axis=1)

    df_validation['val_classindex'] = df_validation.apply(lambda x : match_expectedmatrix_rsv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "index") , axis=1)

    df_validation['val_classcomment'] = df_validation.apply(lambda x : match_expectedmatrix_rsv(x.val_clade, x.val_Gclade, list_expectedmatrix_header, "comment") , axis=1)

    df_validation['val_insertions'] = df_validation.apply(lambda x : seek_insertion_rsv(str(x.nextclade_insertions)) , axis=1)

    df_validation['val_expectedprofile'] = df_validation.apply(lambda x : seek_profile(list_expectedmatrix_change, list_expectedmatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_likelyprofile'] = df_validation.apply(lambda x : seek_profile(list_likelymatrix_change, list_likelymatrix[int(x.val_classindex)][1:]) , axis=1)

    df_validation['val_nocovsub'] = df_validation.apply(lambda x : seek_nocovsub(list_expectedmatrix_change, [y for y in x.stat_missingaapos if y.split(':')[0] in list_target_gene]) , axis=1)

    df_validation['val_nocovins'] = df_validation.apply(lambda x : seek_nocovins(list_expectedmatrix_change, x.stat_missingabspos) , axis=1)

    df_validation['val_missingsub'] = df_validation.apply(lambda x : list_intersect(x.val_expectedprofile, [y for y in (str(x.nextclade_aaSubstitutions) + "," + str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if y not in ['', 'nan']], 'missing', 100) , axis=1)

    df_validation['val_atypicsub'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaSubstitutions)).split(',') if y.split(':')[0] in list_target_gene and y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_atypicindel'] = df_validation.apply(lambda x : list_intersect([y for y in (str(x.nextclade_aaDeletions) + "," + str(x.val_insertions)).split(',') if (y.split(':')[0] in list_target_gene or y[0:3] == 'Ins') and y not in ['', 'nan']], x.val_likelyprofile, 'missing', 100) , axis=1)

    df_validation['val_avisbio'] = df_validation.apply(lambda x : seek_issue_rsv(x.val_glims, x.val_clade, x.summary_qc_seqcontrol, x.val_Fcoverage, x.val_Gcoverage, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_nocovsub, x.val_nocovins, x.val_missingsub, x.nextclade_qc_privateMutations_status, x.summary_vcf_coinf01score, x.summary_vcf_coinf02match, x.summary_vcf_coinf02count, x.summary_vcf_coinf02score, x.val_poi, x.summary_bam_verif, args.coinf_threshold) , axis=1)

    df_validation['val_result'] = df_validation.apply(lambda x : validate_result_rsv(x.val_clade, x.val_dp, x.val_varcount, x.summary_qc_seqcontrol, x.val_avisbio) , axis=1)

    df_validation['val_commentary'] = df_validation.apply(lambda x : gen_commentary_rsv_simple(x.val_classmatch, [y for y in str(x.nextclade_aaSubstitutions).split(',') if y not in ['', 'nan'] and y.split(':')[0] in list_target_gene and y[:-1] in list_profile_base], list_profile_base, [y for y in str(x.nextclade_aaDeletions).split(',') if y not in ['', 'nan'] and y.split(':')[0] in list_target_gene and y[:-1] in list_profile_base], [y for y in str(x.val_insertions).split(',') if y not in ['', 'nan'] and y.rstrip('ACGTN') in list_profile_base], x.val_nocovsub, x.val_nocovins, x.val_classcomment, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_poi), axis=1)

    df_validation['val_commentary_full'] = df_validation.apply(lambda x : gen_commentary_rsv_full(x.val_classmatch, [y for y in str(x.nextclade_aaSubstitutions).split(',') if y not in ['', 'nan'] and y.split(':')[0] in list_target_gene and y[:-1] in list_profile_base], [y for y in str(x.nextclade_aaDeletions).split(',') if y not in ['', 'nan'] and y.split(':')[0] in list_target_gene and y[:-1] in list_profile_base], [y for y in str(x.val_insertions).split(',') if y not in ['', 'nan'] and y.rstrip('ACGTN') in list_profile_base], x.val_nocovsub, x.val_nocovins, x.val_classcomment, x.val_atypicsub.split(';'), x.val_atypicindel.split(';'), x.val_poi), axis=1)

    df_validation['val_error'] = df_validation.apply(lambda x : 'ERROR' if x.val_plate in list_errorplate or x.name in list_errorlist else "OK", axis=1)

    df_validation['stat_simpleprofile'] = df_validation.apply(lambda x : x.val_clade + ' - ' + x.val_atypicsub , axis=1)
    
    df_var_major = filter_dfrownames(df_var, ['|MAJOR'], 'keep')
    df_var_major_indel = filter_dfrownames(df_var_major, ['|DEL|', '|INS|', '|COMPLEX|'], 'keep')

    df_validation_sample = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'drop')
    df_validation_control = filter_dfrownames(df_validation, ['Tpos', 'PCR', 'NT', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')
    df_validation_nt = filter_dfrownames(df_validation, ['NT', 'PCR', 'Neg', 'neg', 'T-', 'TVide', 'Tvide'], 'keep')

    df_var_major_indel[df_var_major_indel > 0].count(axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_indel.tsv', sep='\t', index_label='variant', header=["count"])

    df_validation_control[['summary_consensus_perccoverage', 'val_Gcoverage', 'nextclade_G_clade', 'nextclade_totalAminoacidSubstitutions']].replace(r'^\s*$', np.nan, regex=True).fillna('ININT').to_csv(args.outdir + '/' + run_date + '_' + datatype + '_control.tsv', sep='\t', index_label='control')

    df_validation[['val_result', 'val_commentary', 'val_avisbio', 'val_nocovsub', 'val_nocovins', 'val_missingsub', 'val_atypicsub', 'val_atypicindel', 'summary_bam_meandepth', 'summary_consensus_perccoverage', 'val_Fcoverage', 'val_Gcoverage', 'summary_vcf_dpcount', 'summary_vcf_af05-50count', 'summary_vcf_coinf01match', 'summary_vcf_coinf01count', 'summary_vcf_coinf01score', 'summary_vcf_coinf02match', 'summary_vcf_coinf02count', 'summary_vcf_coinf02score', 'summary_run_id']].replace(r'^\s*$', np.nan, regex=True).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_shortval.tsv', sep='\t', index_label='sample_id')

    df_validation_sample['stat_simpleprofile'].value_counts().sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_runsummary.tsv', sep='\t', index_label='Clade')

    df_metrics = pd.concat([
        df_validation_sample.loc[(df_validation_sample['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') < args.cov_minok) | (df_validation_sample['val_Gcoverage'].apply(pd.to_numeric, errors='coerce') < args.cov_minok)]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_dp'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample.loc[(~df_validation_sample['val_varcount'].isin(['NO','NA'])) & (~df_validation_sample['val_commentary'].isin(['DSC','NA']))]['val_plate'].value_counts(),
        df_validation_sample[df_validation_sample['summary_qc_seqcontrol'].isin(['FAILED'])]['val_plate'].value_counts(),
        df_validation_nt.loc[(df_validation_nt['summary_consensus_perccoverage'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg) | (df_validation_nt['val_Gcoverage'].apply(pd.to_numeric, errors='coerce') >= args.cov_maxneg)]['val_plate'].value_counts()],axis=1).fillna(0).astype(int)
    df_metrics.columns = ['sample_cov<=' + str(args.cov_minok), 'sample_dp>=' + str(args.dp_threshold), 'sample_varcount>=' + str(args.varcount_threshold), 'sample_seq_failed', 'tneg_cov>=' + str(args.cov_maxneg)]
    df_metrics.sort_index(axis = 0).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_metric.tsv', sep='\t', index_label='plate')

    pd.concat([df_validation[df_validation.filter(regex='summary_').columns], df_validation[df_validation.filter(regex='val_').columns], df_validation[df_validation.filter(regex='nextclade_').columns]], axis=1).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_validation.tsv', sep='\t', index_label='sample_id')
    #df_validation_sample_notna = df_validation_sample[~df_validation_sample['summary_consensus_perccoverage_S4'].isin(['NA'])]

    df_export = pd.DataFrame(df_validation_sample, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "VRSMETA"
    df_export['AssayResult'] = df_validation_sample.apply(lambda x : x.val_result.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False)

    df_export = pd.DataFrame(df_validation_sample, columns = ['val_glims'])
    df_export['Instrument ID(s)'] = "NB552333"
    df_export['Analysis authorized by'] = "laurence.josset@chu-lyon.fr"
    df_export['AssayResultTargetCode'] = "COMVRS"
    df_export['AssayResult'] = df_validation_sample.apply(lambda x : x.val_commentary.replace(',',';') if x.val_error != 'ERROR' else 'ERREUR-' + x.val_result.replace(',',';'), axis=1)
    df_export.rename(columns={df_export.columns[0]: 'Sample ID'}).to_csv(args.outdir + '/' + run_date + '_' + datatype + '_fastfinder.csv', sep=',', index = False, mode='a', header=False)
