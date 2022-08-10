import streamlit as st
import pydna
from Bio import Entrez, SeqIO
import requests, sys #
from pydna.dseqrecord import Dseqrecord
from pydna.seqrecord import SeqRecord
from pydna import tm
from pydna.design import primer_design
from Bio.Seq import Seq
from Bio.Restriction import *
#import json

# DMBT1 primers with HindIII and Acc65I (default 55 c and 4 overhang))
# DMBT1 restriction with HindIII and Acc65I

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

def gene_cds(gene,tax):
  gene_term = gene + "[Gene Name] AND " + tax + "[Organism]"
  Entrez.email = ''
  h_search = Entrez.esearch(db="nucleotide", term=gene_term, retmax=1)
  records = Entrez.read(h_search)
  h_search.close()
  identifiers = records['IdList']
  h_fetch = Entrez.efetch(db = 'nucleotide', id =identifiers, rettype = 'gb')
  recs = list(SeqIO.parse(h_fetch,'gb'))
  h_fetch.close()
  handle = Entrez.efetch(db="nucleotide", id=identifiers, rettype="gb", retmode="text")
  a = handle.read()
  for line in a.splitlines():
    if line.find(" CDS ") != -1:
      dots = line.find("..")
      p1 = 2
      while line[dots-p1:dots].isdigit() == True:
        p1 += 1
      nd1 = line[dots-p1+1:dots]
      nd2 = line[dots+2:]
      print(nd1)
      print(nd2)
  # add else
  handle.close()
  s = recs[0].seq
  s = s[int(nd1)-1:int(nd2)+1]
  return s

def prot_str(i,n):
    i_start = i
    # delete human from name (check human appearance in names, if yes check if prot starts from upper case)
    r_data = data[::-1]
    r_start = r_data.index(i)
    r_current = r_start
    r_finish = None
    #st.success(r_data[r_start+1:])
    for i in r_data[r_start+1:]:
        if i[0].isupper() == True:
          r_finish = r_current
          #st.success(i)
        elif i[0].isupper() == False and r_finish != None:
          #st.error(i)
          break
        r_current += 1
    #st.success("r_i and r_x")
    #st.success(r_i)
    #st.success(r_x)
    if i_start == "sequence" or i_start == "primers":
        r_protein = r_data[r_start+1:r_start+1+r_finish+1]
    elif i_start == "coding":
        r_protein = r_data[r_start+1:r_start+1+r_finish]
    #st.success(r_protein)
    prot = r_protein[::-1]
    return prot

def prot_seq(protein,tax_id):
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{protein}) AND (taxonomy_id:{tax_id})&fields=protein_name,gene_names,accession&size=1", headers={"Accept": "text/plain; format=fasta"})
  return r.text

def prot_gene_seq(protein,tax_id,tax):
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{protein}) AND (taxonomy_id:{tax_id})&fields=gene_names", headers={"Accept": "text/plain; format=tsv"})
  entities = r.text
  a_entity = entities.split('\n')
  a_gene = a_entity[1].split(' ')
  gene = a_gene[0]
  s = gene_cds(gene,tax)
  return s

def primers():
    return

st.title('Bioxy')

# Text Input
# save the input text in the variable 'name'
# first argument shows the title of the text input box
# second argument displays a default text inside the text input area
result = ""
task = st.text_input("Please enter the task", result)
 
# display the name when the submit button is clicked
# .title() is used to get the input text string
if(st.button('Submit')):
    #st.clean() find right word
    data = task.split()
    hum_i = try_except('human',data)  # find/create antology of uniprot taxonomies
    seq_i = try_except('sequence',data) 
    prim_i = try_except('primers',data) #add primers -> cloning (look pydna)
    for_i = try_except('for',data)
    if hum_i != -1:
        tax = data[hum_i]
        if tax == 'human':
            tax_id = 9606
    else:
        tax = 'human'
        tax_id = 9606
    if seq_i != -1:
        n,d1,d2 = partial(1)
        if data[seq_i-n].isupper() == True:
            st.success("%s %s sequence" % (tax,data[seq_i-n]))
            gene = data[seq_i-n]
            s = gene_seq(gene,tax)
            if d1 != "" and d2 != "":
                s = s[int(d1):int(d2)+1] #check if -1 is needed after int(d1)
            st.code(s)
        elif data[seq_i-n] == "product":
            n,d1,d2 = partial(2)
            st.success('gene_prot_seq')
        elif data[seq_i-1] == "coding":
            n,d1,d2 = partial(2)
            prot = prot_str('coding',3) 
            if d1 != "" and d2 != "":
              prot = prot[:-1]
            protein = ' '.join(prot)
            #st.success("%s %s coding sequence" % (tax,protein))
            s = prot_gene_seq(protein,tax_id,tax)
            # check gene_cds and dmbt1 aa sequence
            if d1 != "" and d2 != "":
              if d1 == "1":
                 s = s[:(int(d2)*3)]
              else:
                 s = s[((int(d1)-1)*3):(int(d2)*3)]
            st.code(s)
        elif data[seq_i-1] == "plasmid":
            n,d1,d2 = partial(2)
            st.success('plasmid_seq')
        else:
            # work on n
            prot = prot_str('sequence',2) 
            protein = ' '.join(prot)
            # add tax argument as second %s
            st.success("%s %s sequence" % (tax,protein))
            s = prot_seq(protein,tax_id)
            # st.success(s)
            ss = s.splitlines()
            sss = ''.join(ss[1:])
            if d1 != "" and d2 != "":
              ssss = sss[int(d1)-1:int(d2)]
              st.code(ssss)
    elif prim_i != -1: #synchronise
        n,d1,d2 = partial(1)
        if data[seq_i-n].isupper() == True:
            gene = data[seq_i-n] #
            s = gene_seq(gene,tax)
            if d1 != "" and d2 != "":
                s = s[int(d1)-1:int(d2)+1]
            st.success("%s %s sequence" % (tax,data[seq_i-n]))
            st.code(s)
            st.success("%s %s primers" % (tax,data[seq_i-n]))
            # def primers
            dna=Dseqrecord(s)
            ampl = primer_design(dna, target_tm=55.0)
            if prim_i + 1 == for_i:
                st.success("for")
            st.success("forward primer of %s %s" % (tax,data[seq_i-n])) # add temperature
            st.code(ampl.forward_primer.seq)
            st.success("reverse primer of %s %s" % (tax,data[seq_i-n])) # add temperature
            st.code(ampl.reverse_primer.seq)
        else:
            prot = prot_str('primers',1) #
            if d1 != "" and d2 != "": #
                prot = prot[:-1] # 
            protein = ' '.join(prot) #
            s = prot_gene_seq(protein,tax_id,tax)
            if d1 != "" and d2 != "":
              if d1 == "1":
                 s = s[:(int(d2)*3)]
              else:
                 s = s[((int(d1)-1)*3):(int(d2)*3)]
            st.success("%s %s %s-%s coding sequence " % (tax,protein,d1,d2))
            st.code(s)
            st.success("%s %s %s-%s primers" % (tax,protein,d1,d2))
            # def primers
            dna=Dseqrecord(s)
            ampl = primer_design(dna, target_tm=55.0)
            if prim_i + 1 == for_i:
                st.success("for")
            st.success("forward primer of %s %s %s-%s" % (tax,protein,d1,d2)) # add temperature
            st.code(ampl.forward_primer.seq)
            st.success("reverse primer of %s %s %s-%s" % (tax,protein,d1,d2)) # add temperature
            st.code(ampl.reverse_primer.seq)
    else:
        st.error("No result")
