import streamlit as st

st.set_page_config(layout="wide")
st.header("Bioxy is a bioinformatics virtual personal assistant")
st.subheader("with a hands-on experience through a natural language text interface to perform a variety of tasks")
st.write("")
col1, col2, col3 = st.columns(3)

with col1:
    st.image(
            "https://serratus.io/ncbi.png",
            width=200)
with col2:
    st.image("https://upload.wikimedia.org/wikipedia/en/e/e5/UniProt_%28logo%29.png",
            width=200)

with col3:
    st.write("")
    st.write("")
    st.image("https://static.addgene.org/addgene-core/cda5f5fb75/images/common/svg/logo-addgene.svg",
            width=200)
