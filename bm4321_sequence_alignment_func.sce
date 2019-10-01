// BM4321 Genomic Signal Processing
//
// Standalone function file for alignment
//
// Operations
// 1. Global Search
// 2. Local Search
// 3. Semi-Global Search
//
// Objectives
// 1. Familiarization with basic sequence alignment algorithms
//
// Upeka Premaratne, ENTC, University of Moratuwa (upeka@uom.lk) 2016/12/05
// Free to use, distribute and modify for educational purposes with attribution

function g=gap_penalty(k)
    // Gap penalty function
    g_alpha = 0;
    g_beta = 2;
    g = -(g_alpha+g_beta*k);
endfunction

function result=prom_base(base_x,base_y,score_m,score_mm)
   // Match mismatch score
   // score_m - match score
   // score_mm - mismatch score
   if ((base_x==base_y)|(base_x == ascii('A') & base_y == ascii('W'))|(base_x == ascii('T') & base_y == ascii('W'))) then
       result = score_m;
   else
       result = score_mm;
   end
endfunction

function result=comp_base(base_x,base_y,score_m,score_mm)
   // Match mismatch score
   // score_m - match score
   // score_mm - mismatch score
   if ((base_x==base_y)|(base_x == ascii('A') & base_y == ascii('W'))|(base_x == ascii('T') & base_y == ascii('W'))) then
       result = score_m;
   else
       result = score_mm;
   end
endfunction

function matrix_result = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Function to perform global search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //Initialize with gap penalties
    for k_x=1:(len_x+1)
        basic_mat(k_x,1)=gap_penalty(k_x-1);
    end
    
    for k_y=1:(len_y+1)
        basic_mat(1,k_y)=gap_penalty(k_y-1);
    end
    
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]);
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_global(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for global alignment
    score_matrix = scorematrix_global(seq_x,seq_y,score_m,score_mm,gap_penalty);
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    k_y = len_y;
    
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x))]);
        align_y = strcat([pre_y,strrev(ascii(out_y))]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
endfunction

function matrix_result = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y,n] = traceback_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty);
    //disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix); // m - is maximumm value with smallest indices / n - indices [x,y]
    //disp(m);
    //disp(n);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    //disp(k_x);
    //disp(k_y);
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end
    n =[k_x+1,k_y+1];
endfunction

function matrix_result = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+prom_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y 0]); // Add the zero to global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_prom(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty);
    //disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+prom_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

function matrix_result = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    //Function to perform local search
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    
    //No need to initialize because gap penalty will always be negative
    //Recurrance relation
    for k_x=2:(len_x+1)
        for k_y=2:(len_y+1)
            score_match = basic_mat(k_x-1,k_y-1)+comp_base(seq_x(k_x-1),seq_y(k_y-1),score_m,score_mm);
            score_gap_x = basic_mat(k_x,k_y-1)+gap_penalty(1);
            score_gap_y = basic_mat(k_x-1,k_y)+gap_penalty(1);
            basic_mat(k_x,k_y)=max([score_match score_gap_x score_gap_y]); // Same as global search
        end
    end
    matrix_result = basic_mat;
endfunction

function [align_x,align_y] = traceback_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for semi-global alignment
    
    // Assuming x is longer than y to make code simpler
    
    score_matrix = scorematrix_semiglobal(seq_x,seq_y,score_m,score_mm,gap_penalty);
        
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    k_x = len_x;
    [n_max,k_y] = max(score_matrix(len_x+1,:)); 
    k_y=k_y-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+comp_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction

//New modified function
function [align_x,align_y,n] = traceback_prom_modified(seq_x,seq_y,score_m,score_mm,gap_penalty)
    // Traceback for local alignment
    score_matrix = scorematrix_prom(seq_x,seq_y,score_m,score_mm,gap_penalty);
    //disp(score_matrix');
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    
    base_x = seq_x;
    base_y = seq_y;
    
    out_x = [];
    out_y = [];
    
    [m,n]=max(score_matrix);
    k_x = n(1,1)-1;
    k_y = n(1,2)-1;
    
    post_x = ascii(seq_x((k_x+1):len_x));
    post_y = ascii(seq_y((k_y+1):len_y));
    
    // Perform the traceback
    while (k_x>0&k_y>0) then
        k_gx = score_matrix(k_x,k_y+1);
        k_gy = score_matrix(k_x+1,k_y);
        k_diag = score_matrix(k_x,k_y);
        k_c = score_matrix(k_x+1,k_y+1);
        
        if k_c==k_diag+prom_base(seq_x(k_x),seq_y(k_y),score_m,score_mm) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,base_y(k_y)];
            k_x=k_x-1;
            k_y=k_y-1;
        elseif k_c==k_gx+gap_penalty(1) then
            out_x = [out_x,base_x(k_x)];
            out_y = [out_y,45];
            k_x=k_x-1;
        elseif k_c==k_gy+gap_penalty(1) then
            out_x = [out_x,45];
            out_y = [out_y,base_y(k_y)];
            k_y=k_y-1;
        elseif k_c==0;
            break;
        end
    end
    
    // Write the output
    if (k_x>0&k_y==0) then
        pre_x = ascii(base_x(1:k_x));
        pre_y = ascii(45*ones(1,k_x));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);
    elseif (k_x==0&k_y>0) then
        pre_y = ascii(base_y(1:k_y));
        pre_x = ascii(45*ones(1,k_y));
        align_x = strcat([pre_x,strrev(ascii(out_x)),post_x]);
        align_y = strcat([pre_y,strrev(ascii(out_y)),post_y]);       
    else
        align_x = strrev(ascii(out_x));
        align_y = strrev(ascii(out_y));          
    end     
endfunction
