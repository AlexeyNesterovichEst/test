import streamlit as st
#from Bio import Entrez, SeqIO
import Bio
from Bio.Seq import Seq 
seq = Seq("AGCT") 
st.success(seq)
