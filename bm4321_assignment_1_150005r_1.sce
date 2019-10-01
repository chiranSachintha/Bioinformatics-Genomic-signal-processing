// BM4321 Genomic Signal Processing
// Name: C.S Abeygunawardana
// Index No: 150005r
// Assignment 1


clc;
clear all;
close();


exec('bm4321_sequence_alignment_func.sce');
exec('bm4321_gene_prom_region.sce');
exec('bm4321_statistical_alignment_new.sce');

[gp,gn,ncp,ncn]= get_protein_pos_array('ProteinTable.txt');

u_row_p = size(gp,1); 
u_row_n = size(gn,1); 
down_thresh_len = 3; // 3 bases downstream
up_thresh_len = 50;  // 50 bases upstream
safety = 3;

y = ascii('WWWW');
prom_len = 4;

// Question 1
good_p = [];
good_n = [];
located_n = [];
located_p = [];
key_s = [];
key_1 = [];


for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
    
        
        [ax,ay,n] = traceback_local(dna_seq,y,1,-1,gap_penalty);
 
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum

        if x_end < (up_thresh_len-prom_len+1) then
            test_x = dna_seq(x_end+1:x_end+4) //Obtain the promoter
            located_p = [located_p,n_key];
            key_s = [key_s,x_end+1]
            if (((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T')))) then              //Consider only the promoters that match WWWW
                 key_1 = [key_1,x_end+1];
                 good_p = [good_p,n_key];
                 disp(ascii(test_x))

            end
        end
    end
end



for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        disp(ascii(dna_seq));
        
        [ax,ay,n] = traceback_local(dna_seq,y,1,-1,gap_penalty);
        disp(ascii(dna_seq));
        disp(ay);
         
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        if x_end < (up_thresh_len-prom_len+1) then
            test_x = dna_seq(x_end+1:x_end+4) //Obtain the promoter
            located_n = [located_n,n_key];
            key_s = [key_s,x_end+1]
            if (((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T')))) then              //Consider only the promoters that match WWWW
                 key_1 = [key_1,x_end+1];
                 good_n = [good_n,n_key];
                 
            end
        end
    end
end  
 

 
disp("QUESTION 1 answers");

disp(msprintf('Number of Genes of Sense Strand: %d',(u_row_p)))
disp(msprintf('Number of Genes of Anti-Sense Strand: %d',(u_row_n)))
disp(msprintf('Total Number of Genes: %d',(u_row_p+u_row_n)))

disp("LOCAL ALIGNMENT : INTACT QUERY : m=1");
disp(msprintf('Number of Genes with Potential Promoters: %d',length(key_1)))


//disp(key_s);
key_m = key_1;
key_1 = up_thresh_len-key_1+1;
//disp(key_1);
mean_prom_location = (mean(key_1));
disp(mean_prom_location)
promoter_percentage = length(key_1)/(u_row_p+u_row_n)*100;
disp(msprintf('Percentage of Genes with Potential Promoters: %5.3f%%',promoter_percentage))
figure
histplot(0:up_thresh_len, key_1,normalization=%f, style=16);
title('Distribution of Promoters (Intact Query)','fontsize',5)
xlabel('Upstream Position','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3) 
    
// Question 2
key_2 = [];
key_3 = [];
match_len = [];

for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        
        [ax,ay,n] = traceback_local(dna_seq,y,1,-1,gap_penalty);
        
        
        //From scoring matrix
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        
        match = 0
        if x_end < (up_thresh_len-prom_len+1) then
            while (x_end <= (d_len-1)) then //Keep looking for consecutive Ws within the gene starting from the promoter
                if (dna_seq(x_end+1)~= ascii('T')&& dna_seq(x_end+1)~= ascii('A')) then
                    break
                else
                    match = match+1;
                end
                x_end = x_end+1;
            end
            
            if match == length(y) then
                key_2 = [key_2,n(1)];
            end
            if match >= length(y) then
                match_len = [match_len,match];
                key_3 = [key_3,n(1)];
            end
        end
    end
end

for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    //disp(n_key);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        //disp(ascii(dna_seq));
        
        [ax,ay,n] = traceback_local(dna_seq,y,1,-1,gap_penalty);
        
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        
        // Consider only the 50 bases upstream belonging to the coding gene
        match = 0
        if x_end < (up_thresh_len-prom_len+1) then
            while (x_end <= (d_len-1)) then //Keep looking for consecutive Ws within the gene starting from the promoter
                if (dna_seq(x_end+1)~= ascii('T')&& dna_seq(x_end+1)~= ascii('A')) then
                    break
                else
                    match = match+1;
                end
                x_end = x_end+1;
            end
            
            if match == length(y) then
                key_2 = [key_2,n(1)];
            end
            if match >= length(y) then
                match_len = [match_len,match];
                key_3 = [key_3,n(1)];
            end
        end
    end
end

disp("QUESTION 2");
key_2 = up_thresh_len+down_thresh_len-key_2+1;
key_3 = up_thresh_len+down_thresh_len-key_3+1;


figure
histplot(0:20,match_len,normalization=%f, style=16);
title('Consecutive W s Distribution (Intact Query)','fontsize',5)
xlabel('Number of consecutive W s','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3)

disp(msprintf('Maximum number of consecutive Ws (N): %d',max(match_len)))
N = max(match_len);

// Question 3

acc_motif = [];  // Rows of sequences each N bases long
count = 1;
negatives = 0;
row_p = length(good_p); // number of genes positive strand
row_n = length(good_n); // number of genes negative strand

disp(msprintf('From sense strand: %d',row_p));
disp(msprintf('From anti sense strand: %d',row_n));

for n=1:row_p         // consider the sense strand
    n_key = good_p(n);
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    key = key_m(count);
    if (up_thresh_len-key+1) >= N
        dna_row = dna_seq(key:key+N-1);
        d_len = length(dna_row);
        if d_len == N
            save_fasta('promoter.fasta',sprintf('>sequence %i',n_key),dna_row)
            acc_motif = [acc_motif; dna_row];
        end
    else
        negatives = negatives +1
    end
    count = count +1;
end

for n=1:row_n         // consider the anti-sense strand
    n_key = good_n(n);
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    key = key_m(count);
    if (up_thresh_len-key+1) >= N
        dna_row = dna_seq(key:key+N-1);
        d_len = length(dna_row);
        if d_len == N
            disp(ascii(dna_row));
            save_fasta('promoter.fasta',sprintf('>sequence %i',n_key+u_row_p),dna_row)
            acc_motif = [acc_motif; dna_row];
        end
    else
        negatives = negatives +1
    end
    count = count + 1;
end

// Question 3
disp(msprintf('Number of (N) long sequences from Genes Potential Promoter that belong to motif : %d',size(acc_motif,1)))
disp(msprintf('Number of Genes with Potential Promoters that do not belong to motif : %d',negatives))

ppm = get_ppm(acc_motif);
disp("Position Probability Matrix from Intact:");
disp(ppm);

// Question 4

good_p2 = [];
good_n2 = [];
located_n2 = [];
located_p2 = [];
key_s2 = [];
key_1_2 = [];

promoter_lengths = [];


for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    //disp(n_key);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        disp(ascii(dna_seq));
        
        [ax,ay,n] = traceback_local(dna_seq,y,3,-3,gap_penalty);
        disp(ascii(dna_seq));
        disp(ay);
        
        // From scoring matrix
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        if (length(ay)<=d_len) then
            //disp(ascii(dna_seq(x_end+1:length(ay))));
            prom_len2 = length(dna_seq(x_end+1:length(ay)));

            promoter_lengths = [promoter_lengths,prom_len2];
            // Consider only the 50 bases upstream belonging to the coding gene
            if (x_end < (up_thresh_len-prom_len2+1)) then
                //if length(ay)
                test_x = dna_seq(x_end+1:x_end+prom_len2) //Obtain the promoter
                located_p2 = [located_p2,n_key];
                key_s2 = [key_s2,x_end+1]
                if length(test_x)==prom_len then
                    if (((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T')))) then              //Consider only the promoters that match WWWW
                         key_1_2 = [key_1_2,x_end+1];
                         good_p2 = [good_p2,n_key];
                         disp(ascii(test_x))
                     
                     end
                elseif length(test_x)==prom_len+1 then
                    if ((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))) then              //Consider only the promoters that match W-WWW / WW-WW / WWW-W
                         key_1_2 = [key_1_2,x_end+1];
                         good_p2 = [good_p2,n_key];
                     
                     end
                elseif length(test_x)==prom_len+2 then
                    if ((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T'))))) then              //Consider only the promoters that match WWW--W / WW-W-W / WW--WW / W-WW-W / W-W-WW / W--WWW
                         key_1_2 = [key_1_2,x_end+1];
                         good_p2 = [good_p2,n_key];
                         disp(ascii(test_x))
                     
                     end
                end
            
            end
        
        end
    
    end
end


for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        
        [ax,ay,n] = traceback_local(dna_seq,y,3,-3,gap_penalty);

        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        if (length(ay)<=d_len) then
            prom_len2 = length(dna_seq(x_end+1:length(ay)));

            promoter_lengths = [promoter_lengths,prom_len2];
            if (x_end < (up_thresh_len-prom_len2+1)) then
                //if length(ay)
                test_x = dna_seq(x_end+1:x_end+prom_len2) //Obtain the promoter
                located_p2 = [located_p2,n_key];
                key_s2 = [key_s2,x_end+1]
                if length(test_x)==prom_len then
                    if (((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T')))) then              //Consider only the promoters that match WWWW
                         key_1_2 = [key_1_2,x_end+1];
                         good_n2 = [good_n2,n_key];
                       
                     end
                elseif length(test_x)==prom_len+1 then
                    if ((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))) then              //Consider only the promoters that match W-WWW / WW-WW / WWW-W
                         key_1_2 = [key_1_2,x_end+1];
                         good_n2 = [good_n2,n_key];
                        
                     end
                elseif length(test_x)==prom_len+2 then
                    if ((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(2)==ascii('A')))|((test_x(2)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(3)==ascii('A')))|((test_x(3)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T')))))|((((test_x(1)==ascii('A')))|((test_x(1)==ascii('T'))))&(((test_x(4)==ascii('A')))|((test_x(4)==ascii('T'))))&(((test_x(5)==ascii('A')))|((test_x(5)==ascii('T'))))&(((test_x(6)==ascii('A')))|((test_x(6)==ascii('T'))))) then              //Consider only the promoters that match WWW--W / WW-W-W / WW--WW / W-WW-W / W-W-WW / W--WWW
                         key_1_2 = [key_1_2,x_end+1];
                         good_n2 = [good_n2,n_key];
                         disp(ascii(test_x))
                     end
               
                end
           
            end
        
        end
    
    end
end  
 
 


disp(msprintf('Total Number of Genes: %d',(u_row_p+u_row_n)))
disp("LOCAL ALIGNMENT : NON-INTACT QUERY : m=3");
disp(msprintf('Maximum length of promoter in Local Search with Non-Intact-Query: %d',max(promoter_lengths)))
disp(msprintf('Number of Genes with Potential Promoters: %d',length(key_1_2)))

key_m2 = key_1_2;
key_1_2 = up_thresh_len-key_1_2+1;
mean_prom_location2 = (mean(key_1_2));
disp(mean_prom_location2)
promoter_percentage2 = length(key_1_2)/(u_row_p+u_row_n)*100;
disp(msprintf('Percentage of Genes with Potential Promoters: %5.3f%%',promoter_percentage2))
figure
histplot(0:up_thresh_len, key_1_2,normalization=%f, style=16);
title('Distribution of Promoter (Non-Intact Query)','fontsize',4)
xlabel('Upstream Position','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3) 

key_2 = [];
key_3 = [];
key_4 = [];
match_len2 = [];

for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        
        [ax,ay,n] = traceback_local(dna_seq,y,3,-3,gap_penalty);
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        

        match = 0
        if (length(ay)<=d_len) then
            prom_len2 = length(dna_seq(x_end+1:length(ay)));
            // Consider only the 50 bases upstream 
            if (x_end < (up_thresh_len-prom_len2+1)) then
                while (x_end <= (d_len-1)) then //Keep looking for consecutive Ws within the gene starting from the promoter upstream
                    if (dna_seq(x_end+1)~= ascii('T')&& dna_seq(x_end+1)~= ascii('A')) then
                        break
                    else
                        match = match+1;
                    end
                    x_end = x_end+1;
                end
                if match == length(y) then
                    key_2 = [key_2,n(1)];
                end
                if match >= length(y) then
                    match_len2 = [match_len2,match];
                    key_3 = [key_3,n(1)];
                elseif match < length(y) then
                    match_len2 = [match_len2,match];
                    key_3 = [key_3,n(1)];
                end
            end
        end
    end
end

for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        
        [ax,ay,n] = traceback_local(dna_seq,y,3,-3,gap_penalty);

        //From scoring matrix
        y_end = n(2)-1;  //Smallest y index of maximum
        x_end = n(1)-1;  //Smallest x index of maximum
        
        // Consider only the 50 bases upstream belonging to the coding gene
        match = 0
        if x_end < (up_thresh_len-prom_len+1) then
            while (x_end <= (d_len-1)) then //Keep looking for consecutive Ws within the gene starting from the promoter
                if (dna_seq(x_end+1)~= ascii('T')&& dna_seq(x_end+1)~= ascii('A')) then
                    break
                else
                    match = match+1;
                end
                x_end = x_end+1;
            end
            
            if match == length(y) then
                key_2 = [key_2,n(1)];
            end
            if match >= length(y) then
                match_len2 = [match_len2,match];
                key_3 = [key_3,n(1)];
            elseif match < length(y) then
                match_len2 = [match_len2,match];
                key_3 = [key_3,n(1)];
            end
        end
    end
end

key_2 = up_thresh_len+down_thresh_len-key_2+1;
key_3 = up_thresh_len+down_thresh_len-key_3+1;

figure
histplot(0:20,match_len2,normalization=%f, style=16);
title('Consecutive Distribution W s (Non-Intact Query)','fontsize',5)
xlabel('Number of consecutive W s','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3)

disp(msprintf('Maximum number of consecutive Ws (N): %d',max(match_len2)))
N = max(match_len2);

acc_motif2 = [];  // Rows of sequences each N bases long
count2 = 1;
negatives2 = 0;
row_p2 = length(good_p2); // number of genes positive strand
row_n2 = length(good_n2); // number of genes negative strand


disp(msprintf('From sense strand: %d',row_p2));
disp(msprintf('From anti sense strand: %d',row_n2));

for n=1:row_p2         // consider the sense strand
    n_key = good_p2(n);
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    key = key_m2(count2);
    if (up_thresh_len-key+1) >= N
        dna_row = dna_seq(key:key+N-1);
        d_len = length(dna_row);
        if d_len == N
            save_fasta('promoter2.fasta',sprintf('>sequence %i',n_key),dna_row)
            acc_motif2 = [acc_motif2; dna_row];
        end
    else
        negatives2 = negatives2 +1
    end
    count2 = count2 +1;
end

for n=1:row_n2         // consider the anti-sense strand
    n_key = good_n2(n);
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    key = key_m2(count2);
    if (up_thresh_len-key+1) >= N
        dna_row = dna_seq(key:key+N-1);
        d_len = length(dna_row);
        if d_len == N
            save_fasta('promoter2.fasta',sprintf('>sequence %i',n_key+u_row_p),dna_row)
            acc_motif2 = [acc_motif2; dna_row];
        end
    else
        negatives2 = negatives2 +1
    end
    count2 = count2 + 1;
end

disp(msprintf('Number of (N) long sequences from Genes Potential Promoter that belong to motif : %d',size(acc_motif2,1)))
disp(msprintf('Number of Genes with Potential Promoters that do not belong to motif : %d',negatives2))

ppm2 = get_ppm(acc_motif2);
disp("Position Probability Matrix from Non-Intact :");
disp(ppm2);

// Question 5

[w,su]=ppm_info(ppm,[0.25 0.25 0.25 0.25]);
[w2,su2]=ppm_info(ppm2,[0.25 0.25 0.25 0.25]);


fd = mopen('sequence.fasta','rb');
f_out = mopen('results.txt','w'); // Output results
f_out2 = mopen('results2.txt','w'); // Output results
f_out3 = mopen('results3.txt','w'); // Output results

mseek(0,fd);
while (~meof(fd))
    header = mgetl(fd,1);
    keys = strindex(header,' ');
    if (isempty(keys)) then
        break;
    end
    header=ascii(header);
    k_id = ascii(header(2:(keys(1)-1)));
    k_genus = ascii(header((keys(1)+1):(keys(2)-1)));
    k_specie = ascii(header((keys(2)+1):(keys(3)-1)));
    next_entry=0;
    seq = [];
end
mclose(fd);

disp("QUESTION 5");

stat_align_pos = [];
stat_align_score = [];


for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm);
        if (up_thresh_len+1-sm_pos >= N ) then
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            mputl(s_write,f_out);
            stat_align_pos = [stat_align_pos, sm_pos];
            stat_align_score = [stat_align_score, score];
        end
    end
end


for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm);
        if (up_thresh_len+1-sm_pos >= N ) then
            //disp(sm_pos);
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            //disp(s_write);
            mputl(s_write,f_out);
            stat_align_pos = [stat_align_pos, sm_pos];
            stat_align_score = [stat_align_score, score];
        end
    end
end

mclose(f_out);


disp(msprintf('Total Number of Genes: %d',(u_row_p+u_row_n)))

disp("STATISTICAL ALIGNMENT FROM (3)");
disp(msprintf('Number of Genes with Potential Promoters: %d',length(stat_align_pos)))

stat_align_pos = up_thresh_len-stat_align_pos+1;
mean_prom_location3 = (mean(stat_align_pos));
disp(mean_prom_location3)
promoter_percentage3 = length(stat_align_pos)/(u_row_p+u_row_n)*100;
disp(msprintf('Percentage of Genes with Potential Promoters: %5.3f%%',promoter_percentage3))
figure
histplot(0:up_thresh_len,stat_align_pos,normalization=%f, style=16);
title('(a) Distribution of Promoter (Statistical Search (3))','fontsize',5)
xlabel('Upstream Position','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3) 


stat_align_pos2 = [];
stat_align_score2 = [];


for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm2);
        if (up_thresh_len+1-sm_pos >= N ) then
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            mputl(s_write,f_out2);
            stat_align_pos2 = [stat_align_pos2, sm_pos];
            stat_align_score2 = [stat_align_score2, score];
        end
    end
end

for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm2);
        if (up_thresh_len+1-sm_pos >= N ) then
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            mputl(s_write,f_out2);
            stat_align_pos2 = [stat_align_pos2, sm_pos];
            stat_align_score2 = [stat_align_score2, score];
        end
    end
end

mclose(f_out2);

disp("STATISTICAL ALIGNMENT FROM (4)");
disp(msprintf('Number of Genes with Potential Promoters: %d',length(stat_align_pos2)))

stat_align_pos2 = up_thresh_len-stat_align_pos2+1;
mean_prom_location4 = (mean(stat_align_pos2));
disp(mean_prom_location4);
promoter_percentage4 = length(stat_align_pos2)/(u_row_p+u_row_n)*100;
disp(msprintf('Percentage of Genes with Potential Promoters: %5.3f%%',promoter_percentage4))
figure
histplot(0:up_thresh_len,stat_align_pos2,normalization=%f, style=16);
title('(b)Distribution of Promoter(Statistical Search (4))','fontsize',4)
xlabel('Upstream Position','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3) 

stat_align_pos3 = [];
stat_align_score3 = [];

ppm3 = (ppm+ppm2)/2;

disp(ppm3);

//for n=1:row_p2         // consider the sense strand
//    n_key = good_p2(n);
for n_key=1:u_row_p         // consider the sense strand
    dna_seq = get_fasta_at('sequence.fasta',gp(n_key,1)-up_thresh_len,gp(n_key,1)+down_thresh_len+safety,1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm3);
        if (up_thresh_len+1-sm_pos >= N ) then
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            mputl(s_write,f_out3);
            stat_align_pos3 = [stat_align_pos3, sm_pos];
            stat_align_score3 = [stat_align_score3, score];
        end
    end
end

for n_key=1:u_row_n         // consider the anti-sense strand
    dna_seq = get_fasta_at('sequence.fasta',gn(n_key,2)-down_thresh_len+1,gn(n_key,2)+up_thresh_len+safety,-1);
    dna_seq = dna_seq(1:up_thresh_len+down_thresh_len);
    dna_seq = dna_seq(:,$:-1:1);
    d_len = length(dna_seq);
    if (d_len==(up_thresh_len+down_thresh_len)) then
        [s_mat,sm_pos,sm_seq]=get_motif_score(dna_seq,ppm3);
        if (up_thresh_len+1-sm_pos >= N ) then
            //disp(sm_pos);
            if sm_pos~=0 then
                score=log(s_mat(sm_pos));
            else
                score=0;
            end
            s_write = sprintf('%s %s %s %d %2.3f %s %s',k_id,k_genus,k_specie,sm_pos,score,ascii(dna_seq),ascii(sm_seq));
            //disp(s_write);
            mputl(s_write,f_out3);
            stat_align_pos3 = [stat_align_pos3, sm_pos];
            stat_align_score3 = [stat_align_score3, score];
        end
    end
end

mclose(f_out3);

disp("STATISTICAL ALIGNMENT FROM (3)+(4)");
disp(msprintf('Number of Genes with Potential Promoters: %d',length(stat_align_pos3)))

stat_align_pos3 = up_thresh_len-stat_align_pos3+1;
mean_prom_location5 = (mean(stat_align_pos3));
disp(mean_prom_location5)
promoter_percentage5 = length(stat_align_pos3)/(u_row_p+u_row_n)*100;
disp(msprintf('Percentage of Genes with Potential Promoters: %5.3f%%',promoter_percentage5))
figure
histplot(0:up_thresh_len,stat_align_pos3,normalization=%f, style=16);
title('(c)Distribution of Promoter for NZ_CP014692.1(Statistical Search (3)+(4))','fontsize',4)
xlabel('Upstream Position','fontname','times bold','fontsize',3)
ylabel('Genes count','fontname','times bold','fontsize',3) 
