import streamlit as st

# make sidebar
# first
st.title("Documentation")
st.info("[compulsory] & (optional) variables")
option = st.selectbox('Please choose the topic', ('Sequence', 'Primers', 'Restriction', 'Cloning'))
if option == "Sequence":
    c1 = st.container()
    c1.subheader("Gene sequence")
    c1.write("[taxonomy] [GENE NAME] (start-end) sequence")
    c1.markdown("*Example: human DMBT1 1-100 sequence*")

    c2 = st.container()
    c2.subheader("Protein sequence")
    c2.write("[taxonomy] [Protein name] (start-end) sequence")
    c2.markdown("*Example: human Deleted in malignant tumor 1 1-100 sequence*")

    # add partial sequence
    c3 = st.container()
    c3.subheader("Protein coding sequence")
    c3.write("[taxonomy] [Protein name] (start-end) coding sequence")
    c3.markdown("*Example: human Deleted in malignant tumor 1 1-100 coding sequence*")


if option == "Primers":
    c4 = st.container()
    c4.subheader("Gene primers")
    c4.write("[taxonomy] [GENE NAME] (start-end) primers")
    c4.markdown("*Example: human DMBT1 1-100 primers*")
    
    c5 = st.container()
    c5.subheader("Protein coding sequence primers")
    c5.write("[taxonomy] [Protein name] (start-end) primers")
    c5.markdown("*Example: human Deleted in malignant tumor 1 1-100 primers*")
    
if option == "Restriction":
    st.info("In progress")
    
if option == "Cloning":
    st.info("In progress")

