#!/usr/bin/env python 
""" 
This module is for getting summary analyses on adcl, edpl, richness, eveness, Shannon-index alpha diversity and Bray-curtis beta-diversity 
"""

import argparse
from Bio import SeqIO, Entrez, pairwise2
Entrez.email = 'hongyingsun1101@gmail.com'
from Bio.SeqRecord import SeqRecord
import re, time
import os, sys, glob
import pandas as pd
import numpy as np
import random
import uuid
import numpy as np
import pandas as pd
# from skbio.tree import TreeNode
# from skbio import read
# from skbio.stats.distance import DistanceMatrix
import matplotlib as mpl
from scipy import stats
from ast import literal_eval
import sqlite3

""" set up fonts here before importing matplotlib.pylab """
parms = {'font.family': 'serif', 'font.serif': 'Palatino', 'svg.fonttype': 'none'}
mpl.rcParams.update(parms)
from matplotlib import pyplot as pl
import seaborn as sbs
sbs.set(context='talk', style='darkgrid', palette='deep', font='sans-serif')
#pearson R-squared
def pr2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

mock_seqtab=pd.read_csv("CC11.map.SeqTable.csv",index_col=0)
df_sv_list_names=['mock','rdp_10398','rdp_5224','rdp_1017','rdp_92','rdp_12']
taxaFiles=["AbundanceOfTaxIdInSamples_primaryTaxid.csv"]+[d+ "_"+"oneRankEachSV_keepBest.csv" for d in df_sv_list_names[1:]]
taxaDB="taxonomy.db"
conn = create_connection(taxaDB)

nDlists=len(df_sv_list_names)





def percentageCorrect(mock_seqtab, taxaFiles,df_sv_list_names, withMultiplicity=True):
    mock=mock_seqtab.copy()
    #not all tax_ids in the mock are primary
    translate={415850:1463164,195041:45634,592977:1680, 796939:796937, 41791:126333}
    mock.ncbi_tax_id.replace(translate, inplace=True)
    mock.rename(columns={'sourceSeq':'colind'}, inplace=True)
    mock = mock[['community', 'colind','organism', 'ncbi_tax_id', 'multiplicity']]
    mock_gpbySample = mock.groupby('community', as_index=True)
    sample_ids=[list(mock_gpbySample)[i][0] for i in range(len(list(mock_gpbySample)))]
    df_pcorrect=pd.DataFrame(index=sample_ids)
    for i, f in enumerate(taxaFiles[1:],1):
        temp = pd.read_csv(f, index_col=0) #multisample assignment 
        for s in  sample_ids:
            merged=mock_gpbySample.get_group(s).merge(temp[['tax_id', 'tax_name', 'colind']+[s]], how='left', on='colind') #merge on SVs=sourceSeq
            merged.loc[:,"iscorrect"]=(merged['ncbi_tax_id'].astype(str)==merged['tax_id'].astype(str)).astype(int)
            if withMultiplicity:
                df_pcorrect.loc[s,df_sv_list_names[i]]=((merged['iscorrect']*merged['multiplicity'])/merged['multiplicity'].sum()).sum()*100.0 #percentage correct
            else:
                df_pcorrect.loc[s,df_sv_list_names[i]]=(merged['iscorrect']).mean()*100.0 #percentage correct
    return sample_ids, df_pcorrect

sample_ids, df_pcorrect=percentageCorrect(mock_seqtab, taxaFiles,df_sv_list_names)


def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db= "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    
    return (record)
    
    
    
    
# =============================================================================
#     fitched=False
#     while ~fitched:
#         print("hahaha!!!")
#         try:
#             # time.sleep(8)
#             search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
#             fitched=True
#         except:
#             # time.sleep(8)
#             fitched=False
# 
#     return Entrez.read(search)
# =============================================================================

def get_lineage_ids_fromdata(data, uprank):
    """once you have the data from get_tax_data fetch the lineage"""
    #uprank=['kingdom','phylum','class','order','family','genus','species']
    lineage_toparse = data[0]['LineageEx']
#    print("printing data[0]:", data[0]['LineageEx'])
#    print("hahha..............")
    lineage=dict()
    ids=dict()
    for l in lineage_toparse:
#        print("printing l:", l)
        for r in uprank:
            try:
                if l['Rank']==r:
                    lineage[r]=l['ScientificName']
                    ids[r]=l['TaxId']
            except:
                pass
    return lineage, ids
    
def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except ConnectionError as e:
        print(e)

    return None

# this function outputs all the names and ids of every level given a taxid and a connection bulit from SQL database. 
def get_lineage_ids(taxid, conn):
    """ This function gets the names and ids if all parents of the given id """

    query="SELECT nd.tax_id, nd.parent_id, nd.rank, na.tax_id, na.tax_name, na.name_class from nodes nd inner join names na on nd.tax_id=na.tax_id where na.name_class=='scientific name' AND na.tax_id==" + "'"+taxid+"'"

    df = pd.read_sql_query(query, conn)
    #print(df)

    df.columns=['tax_id', 'parent_id', 'rank', 'tax_id_drop', 'tax_name', 'name_class']
    df.drop("tax_id_drop",axis=1, inplace=True)
    rankorder=np.array(['no_rank','superkingdom','phylum','class','order','family','genus','species'])[::-1]
    #print(df['rank'])
    if not df['rank'].iloc[0].strip() in rankorder:
        rankorder=np.append(rankorder,df['rank'].iloc[0])
    rank_ind=np.where(df['rank'].iloc[0]==rankorder)[0][0]
#    if len(df.tax_name.iloc[0].split(" "))>=2:
#        lineage={rankorder[rank_ind]:" ".join(df.tax_name.iloc[0].split(" ")[1:])}
#    else:
#        lineage={rankorder[rank_ind]:df.tax_name.iloc[0]}

    lineage={rankorder[rank_ind]:df.tax_name.iloc[0]}
    ids={rankorder[rank_ind]: taxid}
    stop=False
    temp=df.copy()
    #print(lineage, ids)
    while not stop:
        parent_id=temp['parent_id'].iloc[0]
        if parent_id is None or parent_id=="" or parent_id=='0':
            stop=True
            break
        #print(parent_id)
        query="SELECT nd.tax_id, nd.parent_id, nd.rank, na.tax_id, na.tax_name, na.name_class from nodes nd inner join names na on nd.tax_id=na.tax_id where na.name_class=='scientific name' AND na.tax_id==" + "'"+parent_id+"'"
        temp = pd.read_sql_query(query, conn)
        temp.columns=['tax_id', 'parent_id', 'rank', 'tax_id_drop', 'tax_name', 'name_class']
        temp.drop("tax_id_drop",axis=1, inplace=True)
        lineage.update({temp['rank'].iloc[0]:temp['tax_name'].iloc[0]})

        ids.update({temp['rank'].iloc[0]: temp.tax_id.iloc[0]})



    return lineage, ids

#this function creats a list of files containing mock and the 5 rdp data in the same list with the same format.
def prepareFiles():
    df_sv_list=[]
    for i, f in enumerate(taxaFiles):
        if i==0: #mock
            mock=mock_seqtab.copy()
            sample_ids=mock['community'].unique()
            sample_ids=['CC11CM'+str(i) for i in range(sample_ids.shape[0])]
            #not all tax_ids in the mock are primary
            translate={415850:1463164,195041:45634,592977:1680, 796939:796937, 41791:126333}
            mock.ncbi_tax_id.replace(translate, inplace=True)
            mock.rename(columns={'sourceSeq':'colind','organism':'tax_name','ncbi_tax_id':'tax_id'}, inplace=True)
            sv_ids = mock.colind.unique()
            temp=pd.DataFrame(index=sv_ids,columns=['tax_id']+sample_ids)
            for s in sample_ids:
                mock_s = mock[mock.community==s]
                mock_s.set_index('colind', inplace=True)
                temp.loc[mock_s.index, 'tax_id']=mock_s['tax_id']
                temp.loc[mock_s.index, s]=mock_s['multiplicity']
            temp=temp.fillna(0)
            
        else: #analyzed using dada2/pplacer/RDP
            temp=pd.DataFrame(index=sv_ids,columns=['tax_id']+sample_ids)
            temp1=pd.read_csv(f)
            #drop rank, taxa_name and colind (SV index)
            temp1 = temp1.loc[:,['colind','tax_id']+sample_ids]
            temp1.set_index('colind',inplace=True)
            temp.loc[temp1.index,'tax_id']=temp1['tax_id']
            temp.loc[temp1.index,sample_ids]=temp1[sample_ids]
            temp=temp.fillna(0)

        ##very strange tax_id for the mock is not duplicated but when I set the index as tax_id for the mock the index and the mock becomes duplicated!
        df_sv_list.append(temp)
        #print(temp.head())
    return df_sv_list

# generates the file list which has the mock data and the 5 data generated from pplacer. 
df_sv_list = prepareFiles()

dircs=[i+"/" for i in df_sv_list_names[1:]]

def get_parent(taxid, taxaDB):
    parent_rank=False
    parent_taxid=False
    if taxid=='2':
        return parent_taxid, parent_rank
    #create connection:
    conn = create_connection(taxaDB)
    ranks=np.array(['superkingdom','phylum','subphylum','class','subclass','order','suborder', 'family','genus','species', 'subspecies'])
    try:
        lineage, ids = get_lineage_ids(str(taxid), conn)
    except:
        data_t=get_tax_data(str(taxid))
        lineage, ids = get_lineage_ids_fromdata(data_t, ranks) 

    if type(ids)==dict:
        allparents = ids
    elif type(temp[0])==dict:
        allparents = ids[0]
    else:
        print("Caught exception")
        sys.exit(1)

    rankorder=ranks[::-1]
    #handling when the taxid is itself among parent ids
    for r,i in allparents.items():
        if i==taxid: #and r != "superkingdom":
            taxid_rank=r
            ind_=np.where(rankorder==taxid_rank)[0][0]
            parent_rank=rankorder[ind_+1]


    for i, r in enumerate(rankorder):
        #handling when the taxid is itself among parent ids
        if parent_rank and taxid_rank==r:
            continue

        try:
            parent_taxid=allparents[r]
            parent_rank=r
            break
        except:
            pass
    return parent_taxid, parent_rank



def getUniqueSet(alltaxids):
    out=set(alltaxids[0])
    for l in alltaxids[1:]:
        out=out.union(l)
    return out
allsv=[df_sv_list[i].index.astype(str) for i in range(nDlists)]
allt=[df_sv_list[i].tax_id.astype(str) for i in range(nDlists)]
alls=[df_sv_list[i].columns.astype(str) for i in range(nDlists)]
alltaxa=getUniqueSet(allt)
allsvs=getUniqueSet(allsv)
allsamples=getUniqueSet(alls)

#alln=[get_taxa_adcl(d, prefix).name for d  in dircs]
#allnames=getUniqueSet(alln)
def isSVInsample(svid, sample_id, table):
    """ table has indexes as the unique taxa and rows as samples
        if taxid exist in sample and has abundance>0 return True otherwise return False"""
    assert(table.index.name=="sv_id"),"the table has to have sv_id as its index and the index should be named sv_id"
    if np.any(table.index.isin([svid])):
        if table.loc[svid,sample_id]>0.0: #abundance is not zero
            return True
        else:
            return False
    else:
        return False

def istaxIDEqual(svid, sample_id, table1, table2):
    """this function implicitly assumes that svid exists in both tables for the sample_id and it checks if abundances is>0 in both 
       before comparing their tax_ids"""
    if table1.loc[svid, sample_id]>0 and table2.loc[svid, sample_id]>0:
        tax1 = table1.loc[svid, 'tax_id']
        tax2 = table2.loc[svid, 'tax_id']
    elif table1.loc[svid, sample_id]>0:
        tax1 = table1.loc[svid, 'tax_id']
        tax2=np.nan
    elif table2.loc[svid, sample_id]>0:
        tax1=np.nan
        tax2 = table2.loc[svid, 'tax_id']
    else:
        tax1=np.nan
        tax2=np.nan

    if tax1 == tax2:
        return True, tax1, tax2
    else:
        return False, tax1, tax2



def ranks_off(table1, table2, sv_id, sample_id, taxaDB):
    """A measure that calculates how many ranks is the SV in sample_id in table2 is off from 
       that in table1
       table1 and table2 are pandas dataframes with exactly the same indices and columns 
       columns as samples and indices as SVs"""
    ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    isInTable1=isSVInsample(sv_id, sample_id, table1)
    isInTable2=isSVInsample(sv_id, sample_id, table2)
    #SV has abundance in both
    if isInTable1 and isInTable2:
        isTaxaEq, tax1, tax2 = istaxIDEqual(sv_id, sample_id, table1, table2)
        if isTaxaEq: #SV in both and the corresponding tax_id is equal
            return 0
        else: #SV in both and the corresponding tax_ids differ 
            foundParent=False
            #try parents of tax1
            dummy_id = tax1
            output1=0
            while not foundParent and output1<7:
                print(dummy_id)
                parent_id, parent_rank =get_parent(dummy_id, taxaDB)
                if parent_id:
                    foundParent = parent_id==tax2   #istaxaInsample(parent_id, sample_id, table2)
                    dummy_id = parent_id
                    output1+=1
                else:
                    foundParent=True

            #try parents of tax1
            dummy_id = tax1
            output2=0
            while not foundParent and output2<7:
                print(dummy_id)
                parent_id, parent_rank =get_parent(dummy_id, taxaDB)
                if parent_id:
                    foundParent = parent_id==tax2   #istaxaInsample(parent_id, sample_id, table2)
                    dummy_id = parent_id
                    output2+=1
                else:
                    foundParent=True

            if output1>=output2:
                return output1
            else:
                return output2

    elif isInTable1: #not in table2
        isTaxaEq, tax1, tax2 = istaxIDEqual(sv_id, sample_id, table1, table2)
        parent_id, parent_rank =get_parent(tax1, taxaDB)
        if parent_rank in ranks:
            return ranks.index(parent_rank)+2 #plus 2 because the index starts from 0 and it is the parent
        elif parent_rank == "subspecies":
            return 7
        elif parent_rank == "subphylum":
            return 2
        


    elif isInTable2: #not in table1
        isTaxaEq, tax1, tax2 = istaxIDEqual(sv_id, sample_id, table1, table2)
        parent_id, parent_rank =get_parent(tax2, taxaDB)
        if parent_rank in ranks:
            return ranks.index(parent_rank)+2 #plus 2 because the index starts from 0 and it is the parent
        elif parent_rank == "subspecies":
            return 7
        elif parent_rank == "subphylum":
            return 2
    else: #not in both
        return 0

def get_ranksoff(allsvs, df_sv_list, df_sv_list_names, taxaDB):
    """ a measure that uses table1 taken as the mock and table2 taken as full, setA1, setB1, setC1
        one at a time. The taxa tables table1 and table2 would have the same tax_ids aligned 
        in the index and the same samples aligned in columns """
    mock=df_sv_list[0]
    df_merge=pd.DataFrame(index=list(allsvs))
    df_merge.index.name="sv_id"
    mock_tab1=df_merge.merge(mock,how='left', left_index=True, right_index=True)
    mock_tab1=mock_tab1.fillna(0)
    df_ranksoff_all=pd.DataFrame(index=mock_tab1.columns.values[1:])
    for i, df in enumerate(df_sv_list[1:],1):
        ana_tab2=df_merge.merge(df, how="left", left_index=True, right_index=True)
        ana_tab2=ana_tab2.fillna(0)
        dummy_df=pd.DataFrame(index=mock.index,columns=mock.columns.values)
        for sample_id in mock.columns.values[1:]: #exclude tax_id column
            for sv_id in mock.index.values:
                ranksOff =  ranks_off(mock_tab1, ana_tab2, sv_id, sample_id, taxaDB)
                RAM = mock.loc[sv_id, sample_id]/mock.loc[:, sample_id].sum() #mock relative abundance for taxa in this sample
                dummy_df.loc[sv_id,sample_id]=RAM*ranksOff #ranksOff times relative abundance in mock

            df_ranksoff_all.loc[sample_id,df_sv_list_names[i]]=dummy_df.loc[:,sample_id].sum()
            #df_ranksoff_all.to_csv("ranksOff_bySample.csv")
    fig, ax = pl.subplots(figsize=(12,12))
    sbs.violinplot(data=df_ranksoff_all, ax=ax)
    ax.set_ylabel("Ranks off")
    #fig.savefig("ranksOff_dist_comp.png")

    return df_ranksoff_all

#df_ranksoff_all= get_ranksoff(allsvs, df_sv_list, df_sv_list_names, taxaDB)
    

# generate tables for analysis.
def generateFiles():
    mock=df_sv_list[0]
    df_merge=pd.DataFrame(index=list(allsvs))
    df_merge.index.name="sv_id"
    mock_tab1=df_merge.merge(mock,how='left', left_index=True, right_index=True)
    mock_tab1=mock_tab1.fillna(0)
    ref_tabs={}
    for i, df in enumerate(df_sv_list[1:],1):
        ana_tab2=df_merge.merge(df, how="left", left_index=True, right_index=True)
        ana_tab2=ana_tab2.fillna(0)
        ref_tabs[df_sv_list_names[i]]=ana_tab2
    return ref_tabs

## get the score
#table1=mock_tab1
#table2=ref_tabs["rdp_10398"]
#sample_id="CC11CM0"
#sv_id="AB008213.1.1559"
#
##index error
##    community=CM3
#sv_id="AGXW01000015.688.2204"



def get_score(table1, table2, sv_id, sample_id, taxaDB, rankoff_species_genus = 2, rankoff_genus_family=4, rankoff_family_order=8, rankoff_order_class=16, rankoff_class_phylum=32, accuracyoff=32, species_option=True):
    max_penalty=9999
    rankoff_score=0
    accuracyoff_score=0
    isInTable1 = isSVInsample(sv_id, sample_id, table1)
    isInTable2 = isSVInsample(sv_id, sample_id, table2)
    isTaxaEq, tax1, tax2 = istaxIDEqual(sv_id, sample_id, table1, table2)

    if isInTable1 and isInTable2:
        print("printing tax1:", tax1)
        lineage1, ids1 = get_lineage_ids(str(tax1), conn)
        print("printing tax2:", tax2)
        lineage2, ids2 = get_lineage_ids(str(tax2), conn)
        # table 1 the highest rank is species. 
        if "species" in lineage1.keys(): # in table1,the prediction level is species.
            if "species" in lineage2.keys():
                if lineage1["species"] == lineage1["species"]:
                    rankoff_score +=0
                    accuracyoff_score +=0
                else: # both at species level but different species. 
                    rankoff_score +=0
                    if species_option:
                        accuracyoff_score += accuracyoff
                    else:
                        if lineage1["genus"] == lineage2["genus"]:
                            accuracyoff_score += 0
                        else: # genus are different 
                            accuracyoff_score += accuracyoff
            elif "genus" in lineage2.keys():  # the highest rank in table1 is genus. 
                if species_option:
                    rankoff_score += rankoff_species_genus
                    accuracyoff_score += accuracyoff
                else: # genus option
                    if lineage1["genus"] == lineage2["genus"]:
                        rankoff_score +=0
                        accuracyoff_score +=0
                    else:
                        rankoff_score +=0
                        accuracyoff_score += accuracyoff
            else:
                accuracyoff_score += accuracyoff
                
                if "family" in lineage2.keys():
                    if species_option:
                        offscore = rankoff_species_genus + rankoff_genus_family
                        rankoff_score += offscore
                    else:
                        offscore =  rankoff_genus_family
                        rankoff_score += offscore
                elif "order" in lineage2.keys():
                    if species_option:
                        offscore = rankoff_species_genus + rankoff_genus_family+ rankoff_family_order
                        rankoff_score += offscore
                    else:
                        offscore = rankoff_genus_family + rankoff_family_order
                        rankoff_score += offscore            
                elif "class" in lineage2.keys():
                    if species_option:
                        offscore = rankoff_species_genus + rankoff_genus_family+ rankoff_family_order+ rankoff_order_class
                        rankoff_score += offscore
                    else:
                        offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class
                        rankoff_score += offscore            
                elif "phylum" in lineage2.keys():
                    if species_option:
                        offscore = rankoff_species_genus + rankoff_genus_family+ rankoff_family_order+ rankoff_order_class+rankoff_class_phylum
                        rankoff_score += offscore
                    else:
                        offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class+rankoff_class_phylum
                        rankoff_score += offscore                                         
        # table 1 the highest rank is genus. 
        elif "genus" in lineage1.keys(): # in table1,the prediction level is species.
            if "species" in lineage2.keys():
                if species_option:
                    accuracyoff_score += accuracyoff
                    rankoff_score += rankoff_species_genus
                else:
                    rankoff_score += 0
                    if lineage1["genus"] == lineage2["genus"]:
                        accuracyoff_score +=0
                    else:
                        accuracyoff_score += accuracyoff
            elif "genus" in lineage2.keys():  # the highest rank in table1 is genus. 
                if lineage1["genus"] == lineage2["genus"]:
                    rankoff_score +=0
                    accuracyoff_score +=0
                else:
                    rankoff_score +=0
                    accuracyoff_score += accuracyoff
            else:
                accuracyoff_score += accuracyoff
                
                if "family" in lineage2.keys():
                    offscore =  rankoff_genus_family
                    rankoff_score += offscore
                    
                elif "order" in lineage2.keys():
                    offscore =  rankoff_genus_family+ rankoff_family_order
                    rankoff_score += offscore      
                elif "class" in lineage2.keys():
                    offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class
                    rankoff_score += offscore
          
                elif "phylum" in lineage2.keys():
                    offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class+rankoff_class_phylum
                    rankoff_score += offscore
        # table 1 the highest rank is family. 
        else: # for cases with highest rank of family, order, class, phylum. these cases don't exist in mock data.
            accuracyoff_score=max_penalty;
            rankoff_score=max_penalty;
            
        
#        elif "family" in lineage1.keys(): # in table1,the prediction level is species.
#            if "species" in lineage2.keys():
#                if species_option:
#                    rankoff_score += rankoff_species_genus+rankoff_genus_family
#                else:
#                    rankoff_score += rankoff_genus_family
#                    if lineage1["family"] == lineage2["family"]:
#                        accuracyoff_score +=0
#                    else:
#                        accuracyoff_score += accuracyoff
#            elif "genus" in lineage2.keys():  # the highest rank in table1 is genus. 
#                if lineage1["genus"] == lineage2["genus"]:
#                    rankoff_score +=0
#                    accuracyoff_score +=0
#                else:
#                    rankoff_score +=0
#                    accuracyoff_score += accuracyoff
#            else:
#                accuracyoff_score += accuracyoff
#                
#                if "family" in lineage2.keys():
#                    offscore =  rankoff_genus_family
#                    rankoff_score += offscore
#                elif "order" in lineage2.keys():
#                    offscore =  rankoff_genus_family+ rankoff_family_order
#                    rankoff_score += offscore      
#                elif "class" in lineage2.keys():
#                    offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class
#                    rankoff_score += offscore
#          
#                elif "phylum" in lineage2.keys():
#                    offscore = rankoff_genus_family+ rankoff_family_order+ rankoff_order_class+rankoff_class_phylum
#                    rankoff_score += offscore
    else: # if not in both tables
        rankoff_score=np.nan
        accuracyoff_score=np.nan
        
#    elif isInTable1: # in table1 but not table2
#        off_score=rankoff_species_genus + rankoff_genus_family+ rankoff_family_order+ rankoff_order_class+rankoff_class_phylum 
#        rankoff_score += off_score
#        accuracyoff_score += accuracyoff
#    elif  isInTable2: #in table2 but not table1
#        accuracyoff_score=max_penalty;
#        rankoff_score=max_penalty;
#    else: # not exist in both tables
#        rankoff_score +=max_penalty
#        accuracyoff_score +=max_penalty
        
    score=accuracyoff_score+rankoff_score
    return  rankoff_score, accuracyoff_score, score 


def get_community_score(table1, table2,  taxaDB):
    output_table=table2[::]
    output_table=output_table.drop(columns="tax_id")
    community_list= table1.columns[1:]
    sv_list=table1.index
    for sample_id in community_list:
        print("printing sample_id:", sample_id)
        for sv_id in sv_list:
            print("printing sv_id:", sv_id)
            rankoff_score, accuracyoff_score, score =get_score(table1, table2, sv_id, sample_id, taxaDB)
            output_table.loc[sv_id,sample_id]=score
    return output_table
    

score_table=  get_community_score(table1, table2,  taxaDB)

score_table.to_csv("score_table.csv")
  
    

def main():
    # print( get_tax_data(113287))
    # print(get_lineage_ids_fromdata(get_tax_data(113287),['kingdom','phylum','class','order','family','genus','species'] ))
    # print(create_connection("taxonomy.db"))
#    print(prepareFiles())

    table1=mock_tab1
    table2=ref_tabs["rdp_10398"]
    sample_id="CC11CM0"
    sv_id="AGXW01000015.688.2204"
    score_table=  get_community_score(table1, table2,  taxaDB)

    score_table.to_csv("score_table.csv")




def bray_curtis_distance(refindex, community):
    array1 = df_sv_list[0]["CC11CM"+community]
    array2 = df_sv_list[refindex]["CC11CM"+community]
    return(distance.braycurtis(array1,array2))
def bray_curtis_distance():
    d={1:"rdp10398", 2:"rdp5224",3:"rdp1017",4:"rdp92",5:"rdp12"}
    mock=mock_seqtab.copy()
    sample_ids=mock['community'].unique()
    sample_ids_list=['CC11CM'+str(i) for i in range(sample_ids.shape[0])]
    temp=pd.DataFrame(index=sample_ids_list,columns=["rdp10398", "rdp5224","rdp1017","rdp92","rdp12"])
    for i in range(1,6):
        for j in range(sample_ids.shape[0]):
            bcd=bray_curtis_distance(i, str(j))
            temp.loc["CC11CM"+str(j),d[i]]=bcd
    return (temp)
temp.to_csv("Bray-Curtis-Distance-HS.csv")    



if __name__ == '__main__':
	main()




