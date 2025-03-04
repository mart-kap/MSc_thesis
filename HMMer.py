import os
import csv
import numpy as np
import pyhmmer
import pandas as pd
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
import Bio.Seq
from Bio.Blast import NCBIWWW
from Bio import Entrez
from bs4 import BeautifulSoup
import re
from Bio import SearchIO

Perform_HMM = True
Perform_Clustal1 = True
MakeSelection = True
Curator = True
Perform_Clustal2 = True


#Set folder and read organism taxonomy file
main_dir = "/Users/plant/Desktop/Python/3.7/InterPro/6_ACTIN_MART"
DF = pd.read_csv("/Users/plant/Desktop/Python/3.7/InterPro/0_Database/20241115_DatabaseTaxonomy.csv")

if not os.path.exists(f'{main_dir}/1_Output/0_Temp'):
    os.mkdir(f'{main_dir}/1_Output/0_Temp')
if not os.path.exists(f'{main_dir}/1_Output/0_HMMer'):
    os.mkdir(f'{main_dir}/1_Output/0_HMMer')


if Perform_HMM:
    #Make new CSV file with basic information on length and homologous regions with the consensus
    Alignment_Summary = open(f'{main_dir}/1_Output/1_Alignment_Summary.csv', 'w')
    writer = csv.writer(Alignment_Summary)
    row = 'Taxid', 'ID', 'length','Domains','Domains_Global','Start'
    writer.writerow(row)

    #Read list with species to include (I'd suggest to rewrite it to taxID numbers)
    speciesWishes = pd.read_csv(f"/Users/plant/Desktop/Python/3.7/InterPro/6_ACTIN_MART/Slots_Taxonomy_modifiedlist.csv")
    #speciesWishes = pd.read_csv(f"/Users/plant/Desktop/Python/3.7/InterPro/0_Database/Appendix_DatabaseTaxonomy_JHG.csv")

    Exclude = [item for item in speciesWishes['Taxid'].values.tolist() if item not in DF['Taxid'].values.tolist()]

    print(f"{Exclude} not found")
    species = sorted(list(set(DF.loc[DF['Taxid'].isin(speciesWishes['Taxid'].values.tolist())].Taxid.tolist())))

    #Exclude certain taxIDs (if the common name has multiple taxIDs/subspecies)
    #Remove missing species to prevent errors
    species = [sp for sp in species if sp not in Exclude]

    #CREATE ALIGNMENT FOR HMMer PROFILE using ClustalOmega
    # Input file with FASTAs = input_ACTIN.fasta
    clustalo_cline = ClustalOmegaCommandline(infile=f'{main_dir}/input_ACTIN.fasta',
                                             outfile=f'{main_dir}/1_Output/0_Temp/ProfileAlignment.aln', percentid=True,
                                             distmat_out=f'{main_dir}/1_Output/0_Temp/ProfileDisMat.txt', distmat_full=True,
                                             force=True, outfmt='fasta',
                                             outputorder='input-order')  # , distmat_full_iter=True)
    clustalo_cline()

    #READ ALIGNMENT AND MAKE PROFILE
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(f"{main_dir}/1_Output/0_Temp/ProfileAlignment.aln", digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
    msa.name = b"ACTIN"

    #Build HMM profile and write as file (latter = not required)
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    hmm.consensus
    with open(f"{main_dir}/1_Output/0_Temp/HMM.hmm", "wb") as output_file:
        hmm.write(output_file)

    #Get all proteome files
    path = os.path.abspath('/Users/plant/Desktop/Python/3.7/InterPro/0_Database/1_Proteomes')
    files = [entry.path for entry in os.scandir(path) if entry.is_file()]

    #Write all HMM hits (with certain threshold) in a big fasta file
    fasta = open(f'{main_dir}/1_Output/2_FastaSelection_HMM.fasta', 'w')

    #Loop over all species and run HMM
    for sp_index, sp in enumerate(species):
        file = [item for item in files if f"_{sp}.fasta" in item] #9606 = HUMAN
        file = file[0] #Sometimes 2 files are found, always use first one
        #print(file)

        #if not os.path.exists(f"{main_dir}/1_Output/0_HMMer/{sp}.domtbl"):
        #Actual HMM

        #if not os.path.exists(f"{main_dir}/1_Output/0_HMMer/{sp}.domtbl"):
        pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)#,T=50)
        with pyhmmer.easel.SequenceFile(file, digital=True, alphabet=alphabet) as seq_file:
            hits = pipeline.search_hmm(hmm, seq_file)

        #Write HMM results
        with open(f"{main_dir}/1_Output/0_HMMer/{sp}.domtbl", "wb") as f:
            hits.write(f, format="domains")

        #Read proteome as FASTA
        FASTA = SeqIO.index(file, "fasta")

        #Loop over all HMM hits with score > 70, which is an arbitrary theshold for actin
        #Write ACTIN domains per hit. Don't forget that HMM can write multiple HMM hits per protein e.g. 1-14 and 37-375
        #We extract information
        hitnames = [hit.name.decode("utf-8") for hit in hits if hit.score > 70]
        # else:
        #     input = open(f"{main_dir}/1_Output/0_HMMer/{sp}.domtbl", 'rU')
        #     for qresult in SearchIO.parse(input, 'hmmscan3-domtab'):
        #         hits = qresult.hits
        #         hitnames = [f'{hit.id} {hit.description}' for hit in hits if hit.bitscore > 70]


        for hit in hits:
            ACTIN_domain = []
            ACTIN_domain_global = []
            if hit.name.decode("utf-8") in hitnames:

                #Loop over domains and extract amino acid indeces
                for index, domain in enumerate(hit.domains):
                    ACTIN_domain.append(f'{domain.env_from}-{domain.env_to}')
                    ACTIN_domain_global.extend([domain.env_from,domain.env_to])
                ACTIN_domain = ",".join(ACTIN_domain)

                #Second threshold on length
                alignment_length = max(ACTIN_domain_global) - min(ACTIN_domain_global)
                if (alignment_length > 350) and (alignment_length < 450):
                    sequence = str(FASTA[hit.name.decode("utf-8")].seq).replace("*","")
                    M = [i for i, x in enumerate(sequence) if x == "M"]
                    M = [i - min(ACTIN_domain_global) for i, x in enumerate(sequence) if x == "M" and i - min(ACTIN_domain_global) <= 0]
                    if len(M) != 0:
                        sequence = sequence#[min(ACTIN_domain_global)+M[-1]:]
                    else:
                        M = [0]

                    #if len(sequence) < 600:
                    fasta.write(f'>{sp}_{hit.name.decode("utf-8")}\n')
                    fasta.write(f'{sequence}\n\n')

                    row = sp, hit.name.decode("utf-8"), len(sequence), ACTIN_domain,f'{min(ACTIN_domain_global)}-{max(ACTIN_domain_global)}',min(ACTIN_domain_global)+M[-1],
                    writer.writerow(row)

        print(f'{sp_index+1}/{len(species)} | {DF["Species"].loc[DF["Taxid"] == sp].values.tolist()[0]} ({sp}) | {len(hitnames)} actin-likes')


    fasta.close()

if Perform_Clustal1:
#Above all actin sequences are extracted. Below we align them to remove all ARPs etc.
    print("ClustalO-1")
    clustalo_cline = ClustalOmegaCommandline(infile=f'{main_dir}/1_Output/2_FastaSelection_HMM.fasta',
                                             outfile=f'{main_dir}/1_Output/3_Alignment_Fasta_ClustalW.aln', percentid=True,
                                             distmat_out=f'{main_dir}/1_Output/0_Temp/3_HMMSelection_DisMat.txt', distmat_full=True,
                                             force=True, outfmt='fasta',
                                             outputorder='input-order')  # , distmat_full_iter=True)
    clustalo_cline()

    #Here I use the Distance Matrix of this alignment to remove the ARPs
    #The default distance matrix didn't suit my needs and I rewrote the format.
    #From a 2D matrix to "protein1, protein2, homology"

if MakeSelection:
    DistMat2 = pd.read_table(f"{main_dir}/1_Output/0_Temp/3_HMMSelection_DisMat.txt",delimiter='\s+',skiprows=[0],header=None)
    x = np.asarray(DistMat2)
    i, j = np.triu_indices_from(x, k=1)
    v = x[i, j]
    ijv = np.concatenate((i, j, v)).reshape(3, -1).T
    df_ijv = pd.DataFrame(ijv)
    df_ijv[1] = df_ijv[1]-1
    listx = np.array(x[:,0])
    df_ijv[0] = listx[df_ijv[0].values.tolist()]
    df_ijv[1] = listx[df_ijv[1].values.tolist()]

    #Using Marchantia as a reference, I extracted everythin with a >70% homology
    #You might want extend this with multiple species, as in either Marchantia_ACT1 OR HomoSapiens_ACTB > 70%
    REFs = ["403677_tr|D0NRP9|D0NRP9_PHYIT","3197_tr|A0A1B4XVI7|A0A1B4XVI7_MARPO","2053491_tr|A0A3L6JD80|A0A3L6JD80_THOAR"]
    DF_temp = df_ijv.loc[df_ijv[0].isin(REFs) | df_ijv[1].isin(REFs)]
    DF_temp.to_csv(f"{main_dir}/1_Output/0_Temp/4_DisMat_HMMSelection_DF_all.csv",index=False)
    DF_temp = DF_temp.loc[DF_temp[2] > 50]
    DF_temp.to_csv(f"{main_dir}/1_Output/0_Temp/4_DisMat_HMMSelection_DF.csv",index=False)

if Curator:
    #Read all FASTAs (from the HMMer) and remove all non-actin proteins
    DF_temp = pd.read_csv(f"{main_dir}/1_Output/0_Temp/4_DisMat_HMMSelection_DF.csv")
    FASTA2 = SeqIO.index(f"{main_dir}/1_Output/2_FastaSelection_HMM.fasta", "fasta")
    ACCESSIONS = list(set(DF_temp["0"].values.tolist()+DF_temp["1"].values.tolist()))

    Reduced_list = pd.read_csv("/Users/plant/Desktop/Python/3.7/InterPro/6_ACTIN_MART/1_Output/1_Alignment_Summary.csv")
    Reduced_list['label'] = Reduced_list['Taxid'].astype(str)+"_"+Reduced_list['ID']
    Reduced_list = Reduced_list.loc[Reduced_list['label'].isin(ACCESSIONS)]
    Reduced_list.to_csv(f"{main_dir}/1_Output/4_Alignment_Summary.csv")

    Entrez.email = 'jasper.lamers@hotmail.com'

    #Make new folder
    if not os.path.exists(os.path.join(main_dir,"1_Output","0_Blast")):
        os.mkdir(os.path.join(main_dir,"1_Output","0_Blast"))

    Currated = open(f'{main_dir}/Currated.csv', 'w')
    writer = csv.writer(Currated)
    row = 'ID', 'New', 'Sequence',"DifferentvsInput"
    writer.writerow(row)

    #BLAST MODULE
    for index, item in sorted(enumerate(ACCESSIONS)):
        #try:
        if not (str(FASTA2[item].seq)[-2:] == "CF" and len(str(FASTA2[item].seq)) > 370 and len(str(FASTA2[item].seq)) < 380):
            Blast = False
            if not (os.path.exists(os.path.join(main_dir,"1_Output","0_Blast",f"{item}.xml"))):
                Blast = True
            else:
                result_handle = open(os.path.join(main_dir,"1_Output","0_Blast",f"{item}.xml"), 'r').read()
                XMLfile = BeautifulSoup(result_handle, "xml")
                XMLfile_Hits = XMLfile.find_all('BlastOutput_query-len')

                if len(str(FASTA2[item].seq)) != int(XMLfile_Hits[0].text):
                    Blast = True

            if Blast:
                print("BLAST", item, str(FASTA2[item].seq)[-7:], f'{index}/{len(ACCESSIONS)}')
                taxid = item.split("_")[0]
                sequence = str(FASTA2[item].seq)
                result_handle = NCBIWWW.qblast(program="tblastn", database="core_nt", sequence=sequence, hitlist_size=5,entrez_query=f'txid{taxid}[ORGN]')#, format_type="Text")
                result_handle = result_handle.read()
                result_path = os.path.join(main_dir,"0_Blast")
                f = open(os.path.join(main_dir,"1_Output","0_Blast",f"{item}.xml"), "w")
                f.write(result_handle)
                f.close()
        #except:
         #   print(f"FAILED {item}")

    #READ MODULE
    Count_T = 0
    Count_F = 0

    fasta = open(f'{main_dir}/1_Output/5_Selected_FASTAs2_Curated.fasta', 'w')

    #Loop over FASTAs
    for index, item in enumerate(ACCESSIONS):
        #item = "3702_sp|P93738|ACT9_ARATH"
        #item = "4530_Os12t0163700-01"
        if str(FASTA2[item].seq)[-2:] == "CF" and len(str(FASTA2[item].seq)) > 370 and len(str(FASTA2[item].seq)) < 380:
            Count_T += 1
            fasta.write(f'>{item}\n')
            fasta.write(f'{str(FASTA2[item].seq).replace("*", "")}\n\n')
        else:
            Count_F += 1
            print("READ", item, str(FASTA2[item].seq)[-7:], f'{index}/{len(ACCESSIONS)}')

            with open(os.path.join(main_dir,"1_Output","0_Blast",f"{item}.xml"), 'r') as file:
                result_handle = file.read()

            XMLfile = BeautifulSoup(result_handle, "xml")
            XMLfile_Hits = XMLfile.find_all('Hit')
            XMLfile_TopHits = [hit for hit in XMLfile_Hits if float(hit.Hsp_evalue.text) == 0]
            XMLfile_TopHits = [hit for hit in XMLfile_TopHits if int(hit.Hsp_identity.text)+int(hit.Hsp_gaps.text) >= 370 or int(hit.Hsp_identity.text)+int(hit.Hsp_gaps.text) >= int(hit.select('Hsp_align-len')[0].text)-5]

            NCBI2_hit = ""
            if len(XMLfile_TopHits) > 0:
                for hit in XMLfile_TopHits:
                    ID = hit.Hit_id.text
                    ID = ID.split("|")[1]
                    NCBI = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", retmode="xml")
                    NCBI_XML = BeautifulSoup(NCBI.read(), "xml")

                    NCBI2 = NCBI_XML.find_all('GBQualifier')
                    NCBI2_hit = [item for item in NCBI2 if "translation" in str(item)]
                    if len(NCBI2_hit) != 0:
                        NCBI2_hit = str(NCBI2_hit[0].GBQualifier_value.text)

                        different = (str(FASTA2[item].seq) != NCBI2_hit)
                        row = item, ID, NCBI2_hit, different
                        writer.writerow(row)
                        if different:
                            fasta.write(f'>{item}_Curated\n')
                            fasta.write(f'{str(NCBI2_hit).replace("*", "")}\n\n')
                            break

                    else:
                        handle = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta")
                        handle = [r for r in SeqIO.parse(handle, "fasta")][0]

                        T1 = Bio.Seq.translate(handle.seq)
                        T2 = Bio.Seq.translate(handle.seq[1:])
                        T3 = Bio.Seq.translate(handle.seq[2:])

                        T1 = re.findall(r"(M.*?)\*", str(T1))
                        T2 = re.findall(r"(M.*?)\*", str(T2))
                        T3 = re.findall(r"(M.*?)\*", str(T3))
                        comb = T1 + T2 + T3
                        NCBI2_hit = max(comb, key=len)

                        different = (str(FASTA2[item].seq) != NCBI2_hit)
                        row = item, ID, NCBI2_hit, different
                        writer.writerow(row)

                        if different:
                            fasta.write(f'>{item}_Curated\n')
                            fasta.write(f'{str(NCBI2_hit).replace("*", "")}\n\n')
                            break

            if FASTA2[item].seq == NCBI2_hit:
                if len(str(FASTA2[item].seq)) < 500:
                    fasta.write(f'>{item}_|FLAG|\n')
                    fasta.write(f'{str(NCBI2_hit).replace("*", "")}\n\n')
                else:
                    print(item)

    fasta.close()

if Perform_Clustal2:
#Above all actin sequences are extracted. Below we align them to make the tree (in R) and other analyses
    print("ClustalO-2")
    clustalo_cline2 = ClustalOmegaCommandline(infile=f'{main_dir}/1_Output/5_Selected_FASTAs2_Curated.fasta',
                                             outfile=f'{main_dir}/1_Output/6_Alignment_Fasta_ClustalW.aln', percentid=True,
                                             distmat_out=f'{main_dir}/1_Output/0_Temp/6_DisMat.txt', distmat_full=True,
                                             force=True, outfmt='fasta',
                                             outputorder='input-order')  # , distmat_full_iter=True)
    clustalo_cline2()
