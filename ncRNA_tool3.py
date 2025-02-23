import streamlit as st
from Bio import Entrez, SeqIO
from io import StringIO

# Function to fetch metagenomic sequences
def fetch_metagenomic_sequences(query="metagenome", db="nucleotide", retmax=5):
    """
    Fetch metagenomic sequences from NCBI Nucleotide database.

    Parameters:
    - query: Search term for NCBI
    - db: NCBI database to search (default is 'nucleotide')
    - retmax: Number of sequences to fetch

    Returns:
    - List of SeqRecord objects
    """
    Entrez.email = "abhisheknagargoje842@gmail.com"  # Replace with your email
    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]

    if not ids:
        st.warning("No sequences found for the query.")
        return []

    fasta_data = Entrez.efetch(db=db, id=ids, rettype="fasta", retmode="text")
    sequences = list(SeqIO.parse(fasta_data, "fasta"))

    st.success(f"Fetched {len(sequences)} sequences from NCBI.")
    return sequences

# Function to identify ncRNA
def identify_ncRNA(sequences, known_motifs):
    """
    Identify ncRNA sequences by scanning for known motifs.

    Parameters:
    - sequences: List of SeqRecord objects
    - known_motifs: Dictionary of ncRNA types and their associated motifs

    Returns:
    - Dictionary with sequence IDs and identified ncRNA types
    """
    detected_ncRNA = {}
    for seq in sequences:
        for rna_type, motif in known_motifs.items():
            if motif in str(seq.seq):
                detected_ncRNA[seq.id] = rna_type

    return detected_ncRNA

# Streamlit UI
st.title("Metagenomic ncRNA Annotation Tool")
st.sidebar.header("Options")
query = st.sidebar.text_input("Enter NCBI Search Query", "metagenome")
fetch_button = st.sidebar.button("Fetch Sequences")

# Known ncRNA motifs
known_motifs = {"miRNA": "AGCUU", "snoRNA": "CUGA", "lncRNA": "AUGCU"}

if fetch_button:
    sequences = fetch_metagenomic_sequences(query=query)
    if sequences:
        detected_ncRNA = identify_ncRNA(sequences, known_motifs)

        st.subheader("Identified ncRNA Sequences")
        if detected_ncRNA:
            for seq_id, rna_type in detected_ncRNA.items():
                st.write(f"**{seq_id}**: {rna_type}")
        else:
            st.write("No ncRNA motifs detected.")

# File upload feature for user-provided sequences
uploaded_file = st.file_uploader("Upload a FASTA file for annotation", type=["fasta"])

if uploaded_file is not None:
    # Read the file as a text stream
    uploaded_text = uploaded_file.getvalue().decode("utf-8")
    fasta_io = StringIO(uploaded_text)  # Convert text into StringIO for Biopython
    sequences = list(SeqIO.parse(fasta_io, "fasta"))
    detected_ncRNA = identify_ncRNA(sequences, known_motifs)

    st.subheader("Annotated Uploaded Sequences")
    if detected_ncRNA:
        for seq_id, rna_type in detected_ncRNA.items():
            st.write(f"**{seq_id}**: {rna_type}")
    else:
        st.write("No ncRNA motifs detected in uploaded sequences.")
