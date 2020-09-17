        #alignment_file_emboss = st.sidebar.file_uploader("Buraya EMBOSS dosyanızı giriniz", type=["emboss"])
        #alignment_file_clustal = st.sidebar.file_uploader("Buraya CLUSTAL dosyanızı giriniz", type=["clustal"])
        #alignment_file_stockholm = st.sidebar.file_uploader("Buraya STOCKHOLM dosyanızı giriniz", type=["fstockholm"])

from Bio.Align import MultipleSeqAlignment


        #if seq1 and seq2 is not None:
            #alignments = AlignIO.parse([seq1,seq2], "fasta")
            #st.write(alignments)




        #if seq1 and seq2 is not None:

            #sekans1 = AlignIO.read(seq1,"fasta")
            #sekans2 = AlignIO.read(seq2,"fasta")
            #alignment = MultipleSeqAlignment(
            #    [
            #        SeqRecord(Seq(str(sekans1)), "1.Sekans"),
            #        SeqRecord(Seq(str(sekans2)), "2.Sekans")
            #    ]
            #)
            #st.write(alignment)
            #substitution = alignment.substitutions
            #st.text(substitution)






            aligner = Align.PairwiseAligner()
            #alignments = aligner.align("GAACT", "GAT")
            #for a in alignments:
            #    st.text(a)
            aligner.mode = 'local'
            seq1 = "AGAACTC"
            seq2 = "GAACT"
            score = aligner.score(seq1, seq2)
            st.write(score)
            alignments = aligner.align(seq1,seq2)
            for a in alignments:
                st.text(a)

