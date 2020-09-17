import streamlit as st
from Bio.Seq import Seq 
from Bio import SeqIO
from Bio.SeqUtils import GC
import neatbio.sequtils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Data import CodonTable
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
from Bio import pairwise2




def gc_content(seq):
    result = GC(seq)
    return result

def at_content(seq):
    result = float(str(seq).count("A") + str(seq).count("T")) / len(seq) * 100
    return result

def plot(seq, adenin_color, timin_color, guanin_color, sitozin_color):
    nucleotides = ["Adenin", "Timin", "Guanin", "Sitozin"]
    df = pd.DataFrame(seq, index=nucleotides)
    barlist = plt.bar(nucleotides, seq, label="Nucleotides", width=.2)
    barlist[0].set_color(adenin_color)
    barlist[1].set_color(timin_color)
    barlist[2].set_color(guanin_color)
    barlist[3].set_color(sitozin_color)
    
    st.pyplot()

def main():

    st.title("Bioinformatik Uygulaması")

    st.sidebar.header("Lütfen Seçim Yapınız")
    secenekler = ["Giriş","DNA Sekansı", "Çoklu Sekans Hizalama","Hedef Protein Analizi", "Protein Çözünürlüğü"]
    select_box = st.sidebar.selectbox("Yapmak istediğiniz işlem: ", secenekler)

    if select_box == "Giriş":
        pass
    elif select_box == "DNA Sekansı":
        st.subheader("DNA Sekans Analizi")
        seq_dosya = st.sidebar.file_uploader("Lütfen FASTA yapılı dosyanızı girin", type=["FASTA","fa"])

        if seq_dosya is not None:
            dna = SeqIO.read(seq_dosya, "fasta")
            st.write(dna)

            dna_sekansi = dna.seq

            details = st.radio("Detaylar", ("Açıklama", "Sekansı Göster"))
            if details == "Açıklama":
                st.text("")
                st.write(dna.description)
            elif details == "Sekansı Göster":
                st.text("")
                st.write(dna.seq)

            st.text("")
            st.text("")

            st.subheader("Nükleotid bilgisi")

            st.text("")


            if ("M" and "L") in str(dna_sekansi):
                st.write("Nükleotid bilgilerinin hesaplanabilmesi için lütfen bir **DNA sekansı** giriniz!")
                

            else:
                
                adenin = int(str(dna_sekansi).count("A"))
                guanin = int(str(dna_sekansi).count("G"))
                citosin = int(str(dna_sekansi).count("C"))
                timin = int(str(dna_sekansi).count("T"))
                st.write("**Adenin** sayısı = {0} ".format(adenin))
                st.write("**Timin** sayısı = {0} ".format(timin))
                st.write("**Guanin** sayısı = {0} ".format(guanin))
                st.write("**Sitozin** sayısı = {0} ".format(citosin))

                st.text("")
                st.text("")

                if st.checkbox("Grafik üzerinde göster"):
                    adenin_color = st.beta_color_picker('Adenin için renk seçin', "#F50000")
                    timin_color = st.beta_color_picker('Timin için renk seçin', "#00DE2D")
                    guanin_color = st.beta_color_picker('Guanin için renk seçin', "#1A00FF")
                    sitozin_color = st.beta_color_picker('Sitozin için renk seçin', "#000000")
                    
                    numbers = [adenin, timin, guanin, citosin]
                    plot(numbers, adenin_color, timin_color, guanin_color, sitozin_color)
                
                st.text("")

                st.subheader("İçerik oranları")
                st.text("")

                gc_orani = round(gc_content(dna_sekansi), 2)
                at_orani = round(at_content(dna_sekansi),2)

                st.write("**GC** oranı = % {0}".format(gc_orani))
                st.write("**AT** oranı = % {0}".format(at_orani))

                st.text("")

                st.subheader("Protein Sentezi")
                aa = dna_sekansi.translate()
                aa_frekansi = Counter(str(aa))
                st.text("")

                if st.button("Translasyon Tablosunu Görmek için Tıklayınız"):
                    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
                    st.text(standard_table)
                
                st.text("")

                if st.checkbox("Transkripsiyon"):
                    transkribe = dna_sekansi.transcribe()
                    st.write(transkribe)
                elif st.checkbox("Translasyon"):
                    transle = dna_sekansi.translate()
                    st.write(transle)
                elif st.checkbox("Complement"):
                    st.write(dna_sekansi.complement())
                elif st.checkbox("Amino Asit Sıklığı"):
                    st.write(aa_frekansi)
                elif st.checkbox("Amino Asit Grafiği"):
                    st.text("")
                    aa_color = st.beta_color_picker("Buradan Renk Değiştirebilirsin")
                    plt.bar(aa_frekansi.keys(), aa_frekansi.values(), color = aa_color)
                    st.pyplot()
                elif st.checkbox("Tam Amino Asit İsimleri"):
                    st.text("")
                    aa_ismi = str(aa).replace("*", "")
                    aa3 = utils.convert_1to3(aa_ismi)
                    st.write("**Harf Gösterimi**")
                    st.text(aa_ismi)
                    st.write("**************************************")
                    st.write("**Kısaltma Gösterimi**")
                    st.text(aa3)
                    st.write("**************************************")
                    st.write("**Açık İsim Gösterimi**")
                    st.text(utils.get_acid_name(aa3))              

    elif select_box == "Çoklu Sekans Hizalama":
        seq1 = st.sidebar.file_uploader("1.FASTA dosyanızı giriniz", type=["fasta", "fa"])

        seq2 = st.sidebar.file_uploader("2.FASTA dosyanızı giriniz", type=["fasta", "fa"])

        if seq1 and seq2 is not None:
            sekans1 = SeqIO.read(seq1, "fasta")
            sekans2 = SeqIO.read(seq2, "fasta")

            st.text("")
            st.write("**Hizalanmış Durum:**")

            alignments = pairwise2.align.globalxx(sekans1.seq, sekans2.seq)
            st.text(pairwise2.format_alignment(*alignments[0]))


    elif select_box == "Hedef Protein Analizi":
        pass
    elif select_box == "Protein Çözünürlüğü":
        pass


if __name__ == "__main__":
    main()