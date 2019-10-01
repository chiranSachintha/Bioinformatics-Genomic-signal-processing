// BM4321 Genomic Signal Processing
//
// Operations
// 1. Position probability matrix
// 2. Statitical sequence alignment
// 3. Calculation of ppm entropy
//
// Objectives
// 1. Introduction to statistical sequence alignment
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2019/04/03
// Free to use, distribute and modify for educational purposes with attribution

function result = remove_eols(text_in)
    // Remove EOLs of fasta file
    keys = find(text_in==10);
    if (isempty(keys)) then
        result = text_in;
    else
    text_out = text_in(1:(keys(1)-1));
    k_n = length(keys)-1;
    for k=1:k_n
        text_i = text_in((keys(k)+1):(keys(k+1)-1));
        text_out = [text_out, text_i];
    end
    result = [text_out,text_in((keys(k_n+1)+1):length(text_in))];
    end
endfunction

function comp=get_comp(base)
    if (base==65) then
        comp = 84;
    elseif (base==67) then
        comp = 71;
    elseif (base==71) then
        comp = 67;
    elseif (base==84) then
        comp = 65;
    end
endfunction
function gen_code=get_fasta_at(file_name,g_pos,g_end,strand)
   //Estimate the necessary overread to compensate for the EOL charactors of FASTA files
   g_len = g_end-g_pos;
   if g_len>0 then
          n_extra = floor(g_len/70); 
          n_offset = floor(g_pos/70);
          file_details = fileinfo(file_name);
          file_len = file_details(1);
          fd = mopen(file_name,'rb');
          mseek(0,fd);
          header = mgetl(fd,1);
          g_start = length(header);
          mseek(g_start+g_pos+n_offset,fd);
          //g_start = length(header)+1;
          //mseek(g_start+g_pos,fd);
          raw_code = mget(g_len+n_extra,'c',fd);
          mclose(fd);
          code_i = remove_eols(raw_code);
          //disp(ascii(code_i));
          if strand==1 then
                gen_code = code_i;
          elseif strand==-1 then
                //get complementary strand
                len = length(code_i);
                if len>0 then
                    code_c = [];
                    for k=1:len
                        code_c = [code_c,get_comp(code_i(k))];
                        gen_code = code_c;
                    end
                      //gen_code = gen_code(:,$:-1:1);
                else
                    gen_code = [];
                end 
          end
   else
       gen_code = [];
   end
endfunction

function result=compare_oligo(oligo_a,oligo_b)
    if (length(oligo_a)~=length(oligo_b)) then
        result = -1;
    else
        s = sum(abs(oligo_a-oligo_b))
        if (s==0) then
            result = 1;
        else
            result = 0;
        end
    end
endfunction

function bk=get_amino_key(amino)
// Keys of all amino acids and stop codons
    if (amino==ascii('A')) then
        bk=1; // L-Alanine (Ala/A)
    elseif (amino==ascii('C')) then
        bk=2; // L-Cysteine (Cys/C)
    elseif (amino==ascii('D')) then
        bk=3; // L-Aspartic acid (Asp/D)
    elseif (amino==ascii('E')) then
        bk=4; // L-Glutamic acid (Glu/E)
    elseif (amino==ascii('F')) then
        bk=5; // L-Phenylalanine (Phe/F)
    elseif (amino==ascii('G')) then
        bk=6; // Glycine (Gly/G)
    elseif (amino==ascii('H')) then
        bk=7; // L-Histidine (His/H)
    elseif (amino==ascii('I')) then
        bk=8; // L-Isoleucine (Ile/I)
    elseif (amino==ascii('K')) then
        bk=9; // L-Lysine (Lys/K)
    elseif (amino==ascii('L')) then
        bk=10; // L-Leucine (Leu/L)
    elseif (amino==ascii('M')) then
        bk=11; // L-Methionine (Met/M)
    elseif (amino==ascii('N')) then
        bk=12; // L-Asparagine (Asn/ N)
    elseif (amino==ascii('P')) then
        bk=13; // L-Proline (Pro/P)
    elseif (amino==ascii('Q')) then
        bk=14; // L-Glutamine (Gln/Q)
    elseif (amino==ascii('R')) then
        bk=15; // L-Arginine (Arg/R)
    elseif (amino==ascii('S')) then
        bk=16; // L-Serine (Ser/S)
    elseif (amino==ascii('T')) then
        bk=17; // L-Threonine (Thr/T)
    elseif (amino==ascii('V')) then
        bk=18; // L-Valine (Val/V)
    elseif (amino==ascii('W')) then
        bk=19; // L-Tryptophan (Trp/W)'
    elseif (amino==ascii('Y')) then
        bk=20; // L-Tyrosine (Tyr/Y)
    elseif (amino==ascii('+')|amino==ascii('U')|amino==ascii('O')) then
        bk=21; //Stop codons
    end
endfunction

function bk=get_base_key(base)
    // Base key as A=1, C=2, G=3, T=4 (alphabetic)
    if (base==65) then
        bk=1;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=3;
    elseif (base==84) then
        bk=4;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=1; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=3; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=4; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=3; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=1; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=4; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=1; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

function bk=get_ct_key(base)
    // Base key according to the codon table T=1, C=2, A=3, G=4
    if (base==65) then
        bk=3;
    elseif (base==67) then
        bk=2;
    elseif (base==71) then
        bk=4;
    elseif (base==84) then
        bk=1;
    elseif (base==ascii('Y')|base==ascii('N')) then
        bk=2;
    elseif (base==ascii('R')|base==ascii('W')) then
        bk=3;
    // For multiple possibilities assign one for consistency
    // Evenly distributed as much as possible A=3,C=3, G=3, T=2
    elseif (base==ascii('R')) then
        bk=3; // Assign A (A or G)    
    elseif (base==ascii('Y')) then
        bk=2; // Assign C (C or T) 
    elseif (base==ascii('S')) then
        bk=4; // Assign G (C or G) 
    elseif (base==ascii('W')) then
        bk=1; // Assign T (A or T) 
    elseif (base==ascii('K')) then
        bk=4; // Assign G (G or T) 
    elseif (base==ascii('M')) then
        bk=3; // Assign A (A or C) 
    elseif (base==ascii('B')) then
        bk=2; // Assign C (C or G or T) 
    elseif (base==ascii('D')) then
        bk=3; // Assign G (A or G or T) 
    elseif (base==ascii('H')) then
        bk=1; // Assign T (A or C or T) 
    elseif (base==ascii('V')) then
        bk=3; // Assign A (A or C or G)
    elseif (base==ascii('N')) then
        bk=2; // Assign C (any base)
    end
endfunction

function base_hist=get_base_hist(seq)
    hist=ones(1,4)/4;
    l_seq = length(seq);
    for pos=1:l_seq
        key=get_base_key(seq(pos));
        hist(key)=hist(key)+1;
    end
    base_hist=hist;
endfunction

function ppm=get_ppm(motif)
    // Get the position probability matrix
    [m,n]=size(motif);
    hist_mat = [];
    for pos=1:n
       base_hist = get_base_hist(motif(1:m,pos));
       hist_mat = [hist_mat,base_hist'];
    end
    ppm=hist_mat/(m+1);
endfunction

function [w,su]=ppm_info(ppm,p0)
    [m,n]=size(ppm);
    p_mat = [];
    for k=1:m
        p_row = ones(1,n)/p0(k);
        p_mat = [p_mat; p_row];
    end
    w = log(ppm.*p_mat)/log(2);
    su = ppm.*w;
endfunction

function [s_mat,sm_pos,sm_seq]=get_motif_score(seq,ppm)
    l_seq = length(seq);
    [m,n]=size(ppm);
    score=[];
    max_score = 0;
    max_pos = 0;
    max_seq = [];
    for pos=1:(l_seq-n-1)
        sub_seq = seq(pos:(pos+n-1));
        ss_score=1;
        for pos_ss=1:n
            bk = get_base_key(sub_seq(pos_ss));
            ppm_s = ppm(bk,pos_ss);
            ss_score=ss_score*ppm_s;
        end
        score=[score,ss_score];
        if (ss_score>max_score) then
            max_score=ss_score;
            max_pos= pos;
            max_seq=sub_seq;
        end
    end
    s_mat = score;
    sm_pos = max_pos;
    sm_seq = max_seq;
endfunction

//AF410170.1      GATTTAAAAAGGAGTCTATGGGTTTTTAAG
//KU159153.1      GACTTTCCGACTAGCCTATGGTTCTTCCAA
//JQ770303.1      GATTTAAAAAGGAGTCTATGGGTTTTTAAG
//JX661928.1      GATTTTCGGGTTATCCTATGGATATTGAAA
//KU170123.1      TTTTCTCGGAAAACCCTATGGTTCTTTATG
//KM000002.1      GATTTTCAGGCCATTCTATGGTTGTTCAAA
//KY556640.1      TATTTTCATAGAACTCGATGGTTCTTCAAA
//KY595107.1      GATTTTCAGGACATCTTAGGTTTGTTCAAG
//JX661944.1      GATTTTCAAACTACCCTACAGTTGTTCAAG
//HQ995676.1      GATTTTCAGTGCATCCTATGGTTCTTCAAA
//KU159152.1      GATTTTCAGAACAACCTATGGTTGTTCAAA
//AF410164.1      GATTTGAAAAAAAGTATCTGGTTTTTGAAG
//KY556633.1      GATTTTCAGGCCATTCCATGGTTGTTAAAG
//GU266612.1      GATTTTAAGGCAAACCTATGTTTGCTCAAG
//KU170122.1      TTTTCTCGGAAAACCCTATGGTTCTTTATG
//KU170121.1      TTTTCTCGGAAAACCCTATGGTTCTTTATG
//KM085583.1      AATTTTCAGAAGACCCTGTGGTCCTTCAAG
//KM000005.1      GATTTTCAGGCCATTCTATGGTTGTTCAAA

acc_motif2=[ascii('GATTTAAAAAGGAGTCTATGGGTTTTTAAG');
ascii('GACTTTCCGACTAGCCTATGGTTCTTCCAA');
ascii('GATTTAAAAAGGAGTCTATGGGTTTTTAAG');
ascii('GATTTTCGGGTTATCCTATGGATATTGAAA');
ascii('TTTTCTCGGAAAACCCTATGGTTCTTTATG');
ascii('GATTTTCAGGCCATTCTATGGTTGTTCAAA');
ascii('TATTTTCATAGAACTCGATGGTTCTTCAAA');
ascii('GATTTTCAGGACATCTTAGGTTTGTTCAAG');
ascii('GATTTTCAAACTACCCTACAGTTGTTCAAG');
ascii('GATTTTCAGTGCATCCTATGGTTCTTCAAA');
ascii('GATTTTCAGAACAACCTATGGTTGTTCAAA');
ascii('GATTTGAAAAAAAGTATCTGGTTTTTGAAG');
ascii('GATTTTCAGGCCATTCCATGGTTGTTAAAG');
ascii('GATTTTAAGGCAAACCTATGTTTGCTCAAG');
ascii('TTTTCTCGGAAAACCCTATGGTTCTTTATG');
ascii('TTTTCTCGGAAAACCCTATGGTTCTTTATG');
ascii('AATTTTCAGAAGACCCTGTGGTCCTTCAAG');
ascii('GATTTTCAGGCCATTCTATGGTTGTTCAAA')];

////////////////Program Start
//cd /home/upeka/Documents/Teaching/BM4320_Genomic_Signal_Processing/Motif/Code/fasta_files

//h = get_base_hist(acc_motif(:,1));
//ppm = get_ppm(acc_motif);
//disp(ppm);

//[w,su]=ppm_info(ppm,[0.25 0.25 0.25 0.25]);
//disp(sum(su(:,1:12),1));
//disp(su);
//disp(su);          

//fd = mopen('matK_gene_select_plants.fasta','rb');
//f_out = mopen('results.txt','w'); // Output results

//mseek(0,fd);
//while (~meof(fd))
//    header = mgetl(fd,1);
//    //disp(header);
//    keys = strindex(header,' ');
//    if (isempty(keys)) then
//        break;
//    end
//    header=ascii(header);
//    //disp(keys);
//    k_id = ascii(header(2:(keys(1)-1)));
//    k_genus = ascii(header((keys(1)+1):(keys(2)-1)));
//    k_specie = ascii(header((keys(2)+1):(keys(3)-1)));
//    next_entry=0;
//    seq = [];
//    disp(k_id);
//    disp(k_genus);
//    disp(k_specie);
//    next_line = mgetl(fd,1);
//    while (next_line~="")
//        next_line = mgetl(fd,1);
//        if (isempty(next_line)) then
//            break;
//        end
        //next_line=ascii(next_line);
        //disp(next_line);
//        seq=[seq,ascii(remove_eols(next_line))];
        //disp(seq);
//        if (meof(fd)) then
//            break;
//        end
//    end
    //disp(seq);
    
//    [s_mat,sm_pos,sm_seq]=get_motif_score(seq,ppm);
//    disp(sm_pos);
//    if sm_pos~=0 then
//        score=log(s_mat(sm_pos));
//    else
//        score=0;
//    end
//    s_write = sprintf('%s %s %s %d %2.3f %s',k_id,k_genus,k_specie,sm_pos,score,ascii(sm_seq));
//    disp(s_write);
//    mputl(s_write,f_out);
//end

//mclose(fd);
//mclose(f_out);
