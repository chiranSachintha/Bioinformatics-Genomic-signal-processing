// BM4321 Genomic Signal Processing
//
// Operations
// 1. Reading of a FASTA file
// 2. Extraction of coding and non-coding DNA
// 3. Reading of a GenBank Protein Table
// 4. Basic analysis of a coding and non-coding region
//
// Objectives
// 1. Familiarization with the coding regions of different organisms (Archaea, Bacteria, Eukaryota etc.)
// 2. Preliminary investigation of gene promoter regions
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2016/12/05
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

function inc=base_inc(base,c_val)
    if (base==65) then
        inc=c_val+[1 0 0 0]';
    elseif (base==67) then
        inc=c_val+[0 1 0 0]';
    elseif (base==71) then
        inc=c_val+[0 0 1 0]';
    elseif (base==84) then
        inc=c_val+[0 0 0 1]';
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

function result=compare_oligo(oligo_a,oligo_b) // Protiomics
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

function [gene_array_p,gene_array_n,noncoding_array_p,noncoding_array_n]=get_protein_pos_array(filename)
    // Get the coding and non-coding DNA positions from a protein table
    fd = mopen(filename,'r');
    data = mgetl(fd,1);
    ga_p = [];
    ga_n = [];
    
    nca_p = [];
    nca_n = [];
    
    nc_prev_p = 0;
    nc_prev_n = 0;
    
    while (~meof(fd)) 
        data = mgetl(fd,1);
        //disp(type(data));
        keys = strindex(data,ascii(9));
        p_data = strsplit(data,keys);
        pg_start = strtod(p_data(3));
        pg_stop = strtod(p_data(4));
        //disp(strcmp(p_data(5),'-'));
        //disp(strcmp(p_data(5),'+'));
        if (~isempty(p_data(5))) then
            if (strcmp(p_data(5),'-')==1) then
                ga_n = [ga_n; pg_start,pg_stop];
                nca_n = [nca_n; nc_prev_n, (pg_start-1)];
                nc_prev_n = pg_stop+1;
            else             
                ga_p = [ga_p; pg_start,pg_stop];
                nca_p = [nca_p; nc_prev_p, (pg_start-1)];
                nc_prev_p = pg_stop+1;
            end
        end
    end
    mclose(fd);
    gene_array_p=ga_p;
    noncoding_array_p=nca_p;
    gene_array_n=ga_n;
    noncoding_array_n=nca_n;
endfunction

function save_fasta(filename,header_line,gen_code)
    // Save results into a FASTA file
    fd = mopen(filename,'a');
    //fd = mopen(filename,'wc');
    mputl('',fd);
    mputl(header_line,fd);
    gen_len = length(gen_code);
    for key=1:gen_len
        mput(gen_code(key),'c');
        if pmodulo(key,70)==0 then
            mputl('',fd);
        end
    end
    mclose(fd)
endfunction

function [pos] = perfect_match(x,y,n)
   // Check if x contains n continuous permutation of elements in y 
   // x = 'ATTCTTTA';
   // y = 'AT';
   // n = 4;
   xnum = ascii(x);
   ynum = ascii(y);
   pos = [];
   for key=1:length(xnum)-n+1
       key
       test = xnum(key:key+n-1);
       testIny = members(test,ynum)
       if sum(testIny)==n then
            pos = key;
            break
       end
   end
endfunction

//Execution starts here
//exec('bm4321_sequence_alignment_func.sce'); // Load sequence alignment functions
//[gp,gn,ncp,ncn]= get_protein_pos_array('ProteinTable36727_318552.txt');
//u_row = size(gp,1);

//thresh_len = 50;

//y = ascii('TATA');

//key_s = [];

//u_row = 20;
//dna_seq = get_fasta_at('Acetobacter_indonesiensis.fasta',4306,4887,1);
//d_len = length(dna_seq);
//disp('4306 - 4887')
//disp(ascii(dna_seq));
    
//u_row = 5; //To limit the number of alignments
//for n_key=1:u_row
//    dna_seq = get_fasta_at('Acetobacter_indonesiensis.fasta',gp(n_key,1)-3,gp(n_key,1)+50,1);
//    //dna_seq = get_fasta_at('Acetobacter_indonesiensis.fasta',gp(n_key,1),gp(n_key,1)+3,1);
//    d_len = length(dna_seq);
//    disp(n_key);
//    disp(ascii(dna_seq));
    //[ax,ay] = traceback_local(dna_seq,y,1,-1,gap_penalty);
        //keys = strindex(ax,ascii(y));
        //disp(ax);
        //disp(ay);
        //if ~isempty(keys) then
            //key_s = [key_s,keys(1)];
        //end
//end
//disp(mean(key_s));
