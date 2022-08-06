import streamlit as st
#import Bio
import numpy as np
import requests
import pydna
from Bio import Entrez, SeqIO, Seq
#from biopython.Seq import Seq 
#import biopython
seq = Seq("AGCT") 
seq1 = str(seq)
st.success("success")
st.success(seq1)
