import streamlit as st
import sys
import os
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
from Bio import pairwise2, Align
from Bio.Align import substitution_matrices
#from PIL import Image
from chembl_webresource_client.new_client import new_client
import base64
from SessionState import SessionState
#Aşağıdaki libraryler herokuya yüklenmeden önce requirements dosyasına eklenecek
#from rdkit import Chem
#from rdkit.Chem import Descriptors, Lipinski
#import seaborn as sns
from numpy.random import seed
from numpy.random import randn
#from scipy.stats import mannwhitneyu
import json
import pickle
import uuid
import re
from functions import download_button
#import PyPDF2



st.set_option('deprecation.showPyplotGlobalUse', False)

sns.set(style='ticks')

sys.path.append('/usr/local/lib/python3.7/site-packages/')

st.set_option('deprecation.showfileUploaderEncoding', False)

session_state = SessionState.get(name="", button_sent=False)


def mannwhitney(descriptor, verbose=False):

    # seed the random number generator
    seed(1)

    # actives and inactives
    selection = [descriptor, 'bioactivity_class']
    df = df_class[selection]
    active = df[df.bioactivity_class == 'active']
    active = active[descriptor]

    selection = [descriptor, 'bioactivity_class']
    df = df_class[selection]
    inactive = df[df.bioactivity_class == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    #print('Statistics=%.3f, p=%.3f' % (stat, p))

    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
    
    results = pd.DataFrame({'Descriptor':descriptor,
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    filename = 'mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)

    return results

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x

def norm_value(input):
    norm = []
    for i in input['standard_value']:
        i = float(i)
        i = int(i)
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
        
    return x

# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation
def lipinski(smiles, verbose=False):


    moldata= []
    for elem in smiles:

        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
      desc_MolWt = Descriptors.MolWt(mol)
      desc_MolLogP = Descriptors.MolLogP(mol)
      desc_NumHDonors = Lipinski.NumHDonors(mol)
      desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
          
      row = np.array([desc_MolWt,
                      desc_MolLogP,
                      desc_NumHDonors,
                      desc_NumHAcceptors])   
  
      if(i==0):
          baseData=row
      else:
          baseData=np.vstack([baseData, row])
      i=i+1      
  
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

    #You can check it out rdkit library to further reading for these functions

def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    pdf = PyPDF2.PdfFileReader(df)
    b64 = base64.b64encode(pdf.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/pdf;base64,{b64}" download="CSV_Dosyası.csv">CSV Dosyasını Bilgisayarına Yükle</a>'

    return href

def get_binary_file_downloader_html(bin_file, file_label='File'):

    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a>'
    return href

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

        st.warning("Lütfen DNA Sekansınızı sol taraftaki barda bulunan dosya yükleme kısmından içeri aktarınız.")

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

        st.warning("Lütfen karşılaştırma yapacağınız sekansları sol taraftaki barda bulunan dosya yükleme kısmından içeri aktarınız.")
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
        try:
            st.text("")
            arama = st.sidebar.text_input("Aramak İstediğiniz Proteini Yazınız.", "coronavirus")
            target = new_client.target
            target_query = target.search(arama) #Your target protein name that you want to search
            targets = pd.DataFrame.from_dict(target_query)
            st.write("**Chembl Verisi**:")
            st.write(targets)
            st.text("")
            
            hedef_protein = st.number_input("İleri araştırma yapmak istediğiniz Single Proteinin solunda yazan numarayı girin.", min_value=0 ,value=4, format="%d")
            selected_target = targets.target_chembl_id[hedef_protein]
            st.text("")
            st.write("Seçtiğiniz proteinin **ChEMBL ID**'si: {0}".format(selected_target))
            activity = new_client.activity
            res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
            df = pd.DataFrame.from_dict(res)
            
            st.text("")
        except:
            st.warning("Girilen Değer Geçersiz")
        

        if hedef_protein is not None:
            

            if st.checkbox("Seçtiğiniz Proteinin IC50 değerlerine göre elde edilen verisini görüntülemek için tıklayınız."):
                if df.empty:
                    st.warning("Lütfen bir Single Protein seçiniz!")
                else:
                        
                    st.write(df)
                    st.text("")
                
                    st.markdown(download_button(df, 'IC50_Veri.csv', 'CSV Dosyasını Bilgisayarına Yükle'), unsafe_allow_html=True)
            
            
            st.text("")
            st.text("")
            st.markdown("<h3 style='text-align: center; color: red;'>Seçilen Protein için molekül aktivitesini hesaplayan ML programını çalıştırmak isterseniz aşağıda ki butona tıklayınız.</h3>", unsafe_allow_html=True)
            st.text("")
            
            df2 = df[df.standard_value.notna()]
            bioactivity_class = []

            mol_cid = []
            canonical_smiles = []
            standard_value = []
            for unit in df2.standard_value:
                if float(unit) >= 10000:
                    bioactivity_class.append("inactive")
                elif float(unit) <= 1000:
                    bioactivity_class.append("active")
                else:
                    bioactivity_class.append("intermediate")

            for i in df2.molecule_chembl_id:

                mol_cid.append(i)
            
            for i in df2.canonical_smiles:

                canonical_smiles.append(i)
                
            for i in df2.standard_value:

                standard_value.append(i)
                    
                    

            data_tuples = list(zip(mol_cid, canonical_smiles, standard_value, bioactivity_class))
            df3 = pd.DataFrame( data_tuples,  columns=['molecule_chembl_id', 'canonical_smiles', 'standard_value','bioactivity_class' ])
            st.text("")
            if df.empty:
                st.warning("Lütfen bir Single Protein seçiniz!")
            else:
                if st.checkbox("Moleküler Aktivite Hesapla"):
                    st.text("")
                    st.text("")
                    st.write(df3)
                    st.text("")
                    st.markdown(download_button(df3, 'Genel_Veri.csv', 'CSV Dosyasını Bilgisayarına Yükle'), unsafe_allow_html=True)

                    st.text("")
                    if st.selectbox("Yalnızca Aktif Olanları Göster",("Aktif","")):
                        active_data = (df3.loc[df3['bioactivity_class'] == "active"])
                        st.write(active_data)
                        st.text("")
                        st.markdown(download_button(active_data, 'Aktif_Veri.csv', 'CSV Dosyasını Bilgisayarına Yükle'), unsafe_allow_html=True)

            
                st.text("")
                st.text("")
                st.markdown("<h3 style='text-align: center; color: red;'>Lipinski Tanımlayıcılarını Hesaplamak için aşağıdaki butona tıklayınız.</h3>", unsafe_allow_html=True)
                st.text("")
                
                button_sent = st.checkbox("Lipinski Tanımlayıcıları")
                
                if button_sent:
                    session_state.button_sent = True

                if session_state.button_sent:
                    st.subheader("Lipinski Verisi:")
                    st.write("**MW** = Moleküler Ağırlık")
                    st.write("**LogP** = Molekül Çözünürlüğü")
                    st.write("**NumHDonors** = Hidrojen Bağı Vericileri")
                    st.write("**NumHAcceptors** = Hidrojen Bağı Alıcıları")
                    exploratory_data = df3
                    df_lipinski = lipinski(exploratory_data.canonical_smiles)
                    #st.write(df_lipinski)
                    df_combined = pd.concat([exploratory_data,df_lipinski], axis=1)
                    st.subheader("Birleştirilmiş Veri:")
                    st.write(df_combined)
                    st.markdown(download_button(df_combined, 'Birleştirilmiş_Veri.csv', 'CSV Dosyasını Bilgisayarına Yükle'), unsafe_allow_html=True)
                    df_norm = norm_value(df_combined)
                    #st.write(df_norm)
                    df_final = pIC50(df_norm)
                    st.subheader("IC50'nin pIC50'ye dönüştürülmüş halindeki veri seti:")
                    st.write(df_final)
                    st.markdown(download_button(df_final, 'pIC50_Verisi.csv', 'CSV Dosyasını Bilgisayarına Yükle'), unsafe_allow_html=True)
                    df_class = df_final[df_final.bioactivity_class != "intermediate"]

                    def mannwhitney(descriptor, verbose=False):

                        # seed the random number generator
                        seed(1)

                        # actives and inactives
                        selection = [descriptor, 'bioactivity_class']
                        df = df_class[selection]
                        active = df[df.bioactivity_class == 'active']
                        active = active[descriptor]

                        selection = [descriptor, 'bioactivity_class']
                        df = df_class[selection]
                        inactive = df[df.bioactivity_class == 'inactive']
                        inactive = inactive[descriptor]

                        # compare samples
                        stat, p = mannwhitneyu(active, inactive)
                        #print('Statistics=%.3f, p=%.3f' % (stat, p))

                        # interpret
                        alpha = 0.05
                        if p > alpha:
                            interpretation = 'Same distribution (fail to reject H0)'
                        else:
                            interpretation = 'Different distribution (reject H0)'
                        
                        results = pd.DataFrame({'Descriptor':descriptor,
                                                'Statistics':stat,
                                                'p':p,
                                                'alpha':alpha,
                                                'Interpretation':interpretation}, index=[0])
                        filename = 'mannwhitneyu_' + descriptor + '.csv'
                        results.to_csv(filename)

                        return results

                    st.text("")
                    st.text("")
                    session_state.grafik = st.checkbox("Aktif/İnaktif Molekül Grafiği")
                    session_state.mw = st.checkbox("Moleküler Ağırlık/Çözünürlük Grafiği")
                    session_state.pic50 = st.checkbox("pIC50/Moleküler Aktiflik Grafiği")
                    session_state.logp = st.checkbox("Çözünürlük/Moleküler Aktiflik Grafiği")
                    session_state.donors = st.checkbox("Hidrojen Bağı Vericiler/Moleküler Aktiflik Grafiği")
                    session_state.acceptors = st.checkbox("Hidrojen Bağı Alıcılar/Moleküler Aktiflik Grafiği")

                    if session_state.grafik:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Aktif/İnaktif Molekül Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.countplot(x='bioactivity_class', data=df_class, edgecolor='black')

                        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
                        plt.ylabel('Frequency', fontsize=14, fontweight='bold')
                        
                        st.pyplot()
                        #st.markdown(get_table_download_link(veri), unsafe_allow_html=True)
                        
                        #Buralara PDF indirici eklenecek

                    if session_state.mw:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Moleküler Ağırlık/Çözünürlük Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))
                        sns.scatterplot(x='MW', y='LogP', data=df_class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)

                        plt.xlabel('MW', fontsize=14, fontweight='bold')
                        plt.ylabel('LogP', fontsize=14, fontweight='bold')
                        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
                        st.pyplot()
                        
                        #Buralara PDF indirici eklenecek
                        st.write("**Mann-Whitney U Test Verisi**:")
                        st.write(mannwhitney("MW"))

                    if session_state.pic50:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**pIC50/Moleküler Aktiflik Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_class)

                        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
                        plt.ylabel('pIC50 value', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Verisi**:")
                        st.write(mannwhitney("pIC50"))
                    
                    if session_state.logp:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Çözünürlük/Moleküler Aktiflik Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_class)

                        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
                        plt.ylabel('LogP', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Verisi**:")
                        st.write(mannwhitney("LogP"))
                    
                    if session_state.donors:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Hidrojen Bağı Vericiler/Moleküler Aktiflik Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_class)

                        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
                        plt.ylabel('NumHDonors', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Verisi**:")
                        st.write(mannwhitney("NumHDonors"))

                    if session_state.acceptors:
                        st.write("**********************************")
                        st.text("")
                        st.subheader("**Hidrojen Bağı Alıcılar/Moleküler Aktiflik Grafiği**")

                        plt.figure(figsize=(5.5, 5.5))

                        sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_class)

                        plt.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
                        plt.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')
                        st.pyplot()
                        #Buralara PDF indirici eklenecek

                        st.write("**Mann-Whitney U Test Verisi**:")
                        st.write(mannwhitney("NumHAcceptors"))


                
            
            

    elif select_box == "Protein Çözünürlüğü":
        pass


                    



if __name__ == "__main__":
    main()