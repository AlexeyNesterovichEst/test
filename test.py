import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq 
seq = Seq("AGCT") 
st.success(seq)
