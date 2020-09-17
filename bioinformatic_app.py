import streamlit as st
from Bio.Seq import Seq 
from Bio import SeqIO
from Bio.SeqUtils import GC
#import neatbio.sequtils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
from Bio.Data import CodonTable
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
from Bio import pairwise2, Align
from Bio.Align import substitution_matrices

print(substitution_matrices.load())



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
    secenekler = ["Giriş","DNA Sekansı", "Çoklu Sekans Hizalama", "Sekans Puanlaması", "Hedef Protein Analizi", "Protein Çözünürlüğü"]
    select_box = st.sidebar.selectbox("Yapmak istediğiniz işlem: ", secenekler)

    if select_box == "Giriş":
        st.header("Burası giriş sayfasıdır")
        st.text("")
        st.write("Buraya site hakkında genel bilgilendirmeler girilecektir.")
        st.write("lorem ipsum dolor.")
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

            st.write("**1.Sekans** = {0}".format(sekans1.description))
            st.text("")
            st.write("**2.Sekans** = {0}".format(sekans2.description))
            st.text("")
            st.write("**Hizalanmış Durum:**")

            alignments = pairwise2.align.globalxx(sekans1.seq, sekans2.seq)
            st.text(pairwise2.format_alignment(*alignments[0]))
            st.text("")

            if st.checkbox("BLOSUM62 Puanlamasını Görmek için tıklayınız"):

                secenek = st.radio("Seçenekler",("Hizalanmış Sekans Gösterimi","Tüm Sekans Gösterimi"))

                if secenek == "Hizalanmış Sekans Gösterimi":

                    blosum62 = substitution_matrices.load("BLOSUM62")
                    alignment_blossum = pairwise2.align.localds(sekans1.seq, sekans2.seq, blosum62, -10, -1)
                    st.text(pairwise2.format_alignment(*alignment_blossum[0]))
                    
                elif secenek == "Tüm Sekans Gösterimi":
                    blosum62_2 = substitution_matrices.load("BLOSUM62")
                    full_alignment_blossum = pairwise2.align.localds(sekans1.seq, sekans2.seq, blosum62_2, -10, -1)
                    st.text(pairwise2.format_alignment(*full_alignment_blossum[0], full_sequences=True))

                  #Bu kısımda geliştirme olarak gap-penalty ve gap-extension gibi değerleri kullanıcının değiştirebileceği gibi ayarlayabiliriz. 
                  #Geliştirme olarak sekansları taşımak yerine yazılabilir hale de getirebilirim!!!!!

    elif select_box == "Sekans Puanlaması":

        
        st.text("")
        st.text("")
        st.text("")
        st.subheader("Kendi Hesaplamanızı Yapın:")
        
        seq1_puan = st.text_area("1.Sekans")
        seq2_puan = st.text_area("2.Sekans")
        substitution_matrices.load()
        option = st.selectbox('Kullanmak İstediğiniz Matrix?',("BENNER22", 'BENNER6', 'BENNER74', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'DAYHOFF', 'FENG', 'GENETIC', 'GONNET1992', 'JOHNSON', 'JONES', 'LEVIN', 'MCLACHLAN', 'MDM78', 'NUC.4.4', 'PAM250', 'PAM30', 'PAM70', 'RAO', 'RISLER', 'SCHNEIDER', 'STR'))
        st.write('Seçtiğiniz Matrix:', option)
        try:        
            aligner = Align.PairwiseAligner()
            if option == "BENNER22":
                matrix = substitution_matrices.load("BENNER22")
                st.text(matrix)
                aligner.substitution_matrix = matrix
                score = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score)) 
            elif option == "BLOSUM62":
                matrix2 = substitution_matrices.load("BLOSUM62")
                st.text(matrix2)
                aligner.substitution_matrix = matrix2
                score2 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score2))
            elif option == "BENNER6":
                matrix3 = substitution_matrices.load("BENNER6")
                st.text(matrix3)
                aligner.substitution_matrix = matrix3
                score3 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score3))
            elif option == "BENNER74":
                matrix4 = substitution_matrices.load("BENNER74")
                st.text(matrix4)
                aligner.substitution_matrix = matrix4
                score4 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score4))
            elif option == "BLOSUM45":
                matrix5 = substitution_matrices.load("BLOSUM45")
                st.text(matrix5)
                aligner.substitution_matrix = matrix5
                score5 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score5))
            elif option == "BLOSUM50":
                matrix6 = substitution_matrices.load("BLOSUM50")
                st.text(matrix6)
                aligner.substitution_matrix = matrix6
                score6 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score6))
            elif option == "BLOSUM80":
                matrix7 = substitution_matrices.load("BLOSUM80")
                st.text(matrix7)
                aligner.substitution_matrix = matrix7
                score7 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score7))
            elif option == "BLOSUM90":
                matrix8 = substitution_matrices.load("BLOSUM90")
                st.text(matrix8)
                aligner.substitution_matrix = matrix8
                score8 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score8))
            elif option == "DAYHOFF":
                matrix9 = substitution_matrices.load("DAYHOFF")
                st.text(matrix9)
                aligner.substitution_matrix = matrix9
                score9 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score9))
            elif option == "FENG":
                matrix10 = substitution_matrices.load("FENG")
                st.text(matrix10)
                aligner.substitution_matrix = matrix10
                score10 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score10))
            elif option == "GENETIC":
                matrix11 = substitution_matrices.load("GENETIC")
                st.text(matrix11)
                aligner.substitution_matrix = matrix11
                score11 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score11))
            elif option == "GONNET1992":
                matrix12 = substitution_matrices.load("GONNET1992")
                st.text(matrix12)
                aligner.substitution_matrix = matrix12
                score12 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score12))
            elif option == "JOHNSON":
                matrix13 = substitution_matrices.load("JOHNSON")
                st.text(matrix13)
                aligner.substitution_matrix = matrix13
                score13 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score13))
            elif option == "JONES":
                matrix14 = substitution_matrices.load("JONES")
                st.text(matrix14)
                aligner.substitution_matrix = matrix14
                score14 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score14))
            elif option == "LEVIN":
                matrix15 = substitution_matrices.load("LEVIN")
                st.text(matrix15)
                aligner.substitution_matrix = matrix15
                score15 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score15))
            elif option == "MCLACHLAN":
                matrix16 = substitution_matrices.load("MCLACHLAN")
                st.text(matrix16)
                aligner.substitution_matrix = matrix16
                score16 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score16))
            elif option == "MDM78":
                matrix17 = substitution_matrices.load("MDM78")
                st.text(matrix17)
                aligner.substitution_matrix = matrix17
                score17 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score17))
            elif option == "NUC.4.4":
                matrix18 = substitution_matrices.load("NUC.4.4")
                st.text(matrix18)
                aligner.substitution_matrix = matrix18
                score18 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score18))
            elif option == "PAM250":
                matrix19 = substitution_matrices.load("PAM250")
                st.text(matrix19)
                aligner.substitution_matrix = matrix19
                score19 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score19))
            elif option == "PAM30":
                matrix20 = substitution_matrices.load("PAM30")
                st.text(matrix20)
                aligner.substitution_matrix = matrix20
                score20 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score20))
            elif option == "PAM70":
                matrix21 = substitution_matrices.load("PAM70")
                st.text(matrix21)
                aligner.substitution_matrix = matrix21
                score21 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score21))
            elif option == "RAO":
                matrix22 = substitution_matrices.load("RAO")
                st.text(matrix22)
                aligner.substitution_matrix = matrix22
                score22 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score22))
            elif option == "RISLER":
                matrix23 = substitution_matrices.load("RISLER")
                st.text(matrix23)
                aligner.substitution_matrix = matrix23
                score23 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score23))
            elif option == "SCHNEIDER":
                matrix24 = substitution_matrices.load("SCHNEIDER")
                st.text(matrix24)
                aligner.substitution_matrix = matrix24
                score24 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score24))
            elif option == "STR":
                matrix25 = substitution_matrices.load("STR")
                st.text(matrix25)
                aligner.substitution_matrix = matrix25
                score25 = aligner.score(seq1_puan, seq2_puan)
                st.write("**Sekans Puanı** = {0} ".format(score25))
            
            #aligner.substitution_matrix = matrix
            #score = aligner.score("ACDQ", "ACDQ")
            #st.write(score)
        except:
            st.text("")

            

                   

                

                    




    elif select_box == "Hedef Protein Analizi":
        pass
    elif select_box == "Protein Çözünürlüğü":
        pass


if __name__ == "__main__":
    main()