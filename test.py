import streamlit as st
import numpy as np
import requests
import pydna
from Bio import Entrez, SeqIO
import sys #
#import json

#  + Human Deleted in malignant brain tumors 1 sequence
# Human Deleted in malignant brain tumors 1 coding sequence (display gene as well)
# DMBT1 primers
# DMBT1 primers with HindIII and Acc65I (default 55 c and 4 overhang))
# DMBT1 cloning to pQM with HindIII and Acc65I

# Then demo complete, expand by remaining in Colab version, create Bioxy graph

WEBSITE_API = "https://rest.uniprot.org"

def get_url(url, **kwargs):
  response = requests.get(url, **kwargs);

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response

def try_except(word,data):
  try:
    i = data.index(word)
  except ValueError:
    i = -1
  return i

def partial(n):
  d1 = ""
  d2 = ""
  if data[seq_i-n].find("-") != -1:
    ss_text = data[seq_i-n].split('-')
    n += 1
    for ss in ss_text:
      if ss.isdigit() == True:
        if d1 == "":
          d1 = ss
        elif d2 == "":
          d2 = ss
  return n,d1,d2

def gene_seq(gene,tax):
  gene_term = gene + "[Gene Name] AND " + tax + "[Organism]"
  Entrez.email = ''
  h_search = Entrez.esearch(db="nucleotide", term=gene_term, retmax=1)
  records = Entrez.read(h_search)
  h_search.close()
  identifiers = records['IdList']
  h_fetch = Entrez.efetch(db = 'nucleotide', id =identifiers, rettype = 'gb')
  recs = list(SeqIO.parse(h_fetch,'gb'))
  h_fetch.close()
  s = recs[0].seq
  return s

def prot_str(i,n):
    # delete human from name (check human appearance in names, if yes check if prot starts from upper case)
    r_data = data[::-1]
    r_i = r_data.index(i)
    if i == 'coding' and n == 3:
      r_i += 1
    #continue here!
    r_n = r_i
    r_x = None
    for i in r_data[r_i+1:]:
        if i[0].isupper() == True:
          r_x = r_n
        elif i[0].isupper() == False and r_x != None:
          break
        r_n += 1
    #while r_data[r_n][0].isupper() == True or r_data[r_n].isdigit():
       #r_n +=1
    #st.success("r_i and r_x")
    #st.success(r_i)
    #st.success(r_x)
    r_protein = r_data[r_i+1:r_x+r_i+2]
    prot = r_protein[::-1]
    return prot

def prot_seq(prot,tax_id):
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{prot}) AND (taxonomy_id:{tax_id})&fields=protein_name,gene_names,accession&size=1", headers={"Accept": "text/plain; format=fasta"})
  return r.text

st.title('Bioxy')

# Text Input
# save the input text in the variable 'name'
# first argument shows the title of the text input box
# second argument displays a default text inside the text input area
result = ""
task = st.text_input("Enter Your Task", result)
 
# display the name when the submit button is clicked
# .title() is used to get the input text string
if(st.button('Submit')):
    data = task.split()
    hum_i = try_except('human',data)  # find/create antology of uniprot taxonomies
    seq_i = try_except('sequence',data) 
    prim_i = try_except('primers',data) #add primers -> cloning (look pydna)
    if hum_i != -1:
        tax = data[hum_i]
        if tax == 'human':
            tax_id = 9606
    else:
        tax = 'human'
    if seq_i != -1:
        n,d1,d2 = partial(1)
        if data[seq_i-n].isupper() == True:
            st.success("%s sequence" % data[seq_i-n])
            gene = data[seq_i-n]
            s = gene_seq(gene,tax)
            if d1 != "" and d2 != "":
                s = s[int(d1):int(d2)+1] #check if -1 is needed after int(d1)
                st.success(s)
            else:
                st.success(s)
        elif data[seq_i-n] == "product":
            n,d1,d2 = partial(2)
            st.success('gene_prot_seq')
        elif data[seq_i-1] == "coding":
            n,d1,d2 = partial(2)
            st.success('prot_gene_seq')
        elif data[seq_i-1] == "plasmid":
            n,d1,d2 = partial(2)
            st.success('plasmid_seq')
        else:
            # work on n
            prot = prot_str('sequence',2) 
            protein = ' '.join(prot)
            # add tax argument as second %s
            st.success("%s sequence" % protein)
            s = prot_seq(protein,tax_id)
            ss = s.splitlines()
            sss = ''.join(ss[1:])
            st.success(sss)
    elif prim_i != -1:
        n,d1,d2 = partial(1)
        if data[seq_i-n].isupper() == True:
            st.success('gene_seq_primers')
        else:
            st.success('prot_seq_primers')
            #prot = prot_str('primers',0) #partial n,d1,d2 = partial(m,"prot")
            #protein = ' '.join(prot)
            #st.success(protein)
    else:
        st.error("No result")
