

clc;
clear all;
close();

exec('bm4321_phylogenetic_trees.sce');

proteinName = "Acetobacter aceti";
numOfgenes = 22;

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


//NZ_CP014692.1 Acetobacter aceti
//NC_013209.1/AP011121.1 Acetobacter pasteurianus
//NZ_CP023657.1/CP023657.1 Acetobacter pomorum
//NZ_CP022699.1/CP022699.1 Acetobacter tropicalis
//NZ_CP014687.1/CP014687.1 Acetobacter persici
//NZ_LN606600.1/LN606600.1 Acetobacter senegalensis
//NZ_CP011120.1/CP011120.1 Acetobacter oryzifermentans
//NZ_CP015164.1 Acetobacter ascendens
//NZ_CP015168.1 Acetobacter ascendens
//NZ_CP021524.1 Acetobacter ascendens
//NZ_CP022374.1/CP022374.1 Acetobacter oryzifermentans
//NZ_AP018515.1/AP018515.1 Acetobacter orientalis
//NZ_CP023189.1/CP023189.1 Acetobacter pomorum
//NC_017100.1/AP011128.1 Acetobacter pasteurianus
//NC_017121.1/AP011135.1 Acetobacter pasteurianus
//NC_017125.1/AP011142.1 Acetobacter pasteurianus
//NC_017146.1/AP011149.1 Acetobacter pasteurianus
//NZ_LN609302.1 Acetobacter ghanensis
//NC_017111.1/AP011156.1 Acetobacter pasteurianus
//NC_017150.1/AP011163.1 Acetobacter pasteurianus
//NC_017108.1/AP011170.1 Acetobacter pasteurianus
//NZ_AP014881.1/AP014881.1 Acetobacter pasteurianus
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
          raw_code = mget(g_len+n_extra,'c',fd);
          mclose(fd);
          code_i = remove_eols(raw_code);
          if strand== 1 then
              gen_code = code_i;
              //gen_code=gen_code(1:343);
          else
              //get complementary strand
              len = length(code_i);
              code_c = [];
              for k=1:len
                  code_c = [code_c,get_comp(code_i(k))];
                  gen_code = code_c;
              end
              gen_code = gen_code(:,$:-1:1);
              if (ascii(gen_code(1))== 'A') then 
                  if (ascii(gen_code(2))== 'A') then
                      gen_code=gen_code(2:length(gen_code));
                  end 
               else
                  gen_code=gen_code(2:length(gen_code));
                  end
                      
              //gen_code=gen_code(1:343);
          end
   else
       gen_code = [];
   end

endfunction

AssignmentFolder = "C:\Users\Chiran\Documents\bioinformatics"
BaseFolderPath = "C:\Users\Chiran\Documents\bioinformatics\protein"
files =  findfiles(BaseFolderPath,'*.txt');
fileNames = gsort(files,'g','i');

gene_pos=zeros(numOfgenes,4);
for i=1:numOfgenes
    fileName=BaseFolderPath+"\"+fileNames(i);
    f=mopen(fileName,'r');
    disp(fileName);
    while(meof(f)==0)
       line=mgetl(f,1);
       indx=strindex(line,ascii(9));                            
       line=strsplit(line,indx);
       line=stripblanks(line,%t); //strip leading and trailing blanks including tabs
//           if (line(size(line)(1))==proteinName)||(line(size(line)(1))=='MULTISPECIES: '+proteinName) then
       if (grep(line(size(line)(1)),proteinName)==1)
           gPos(i,1)=strtod(line(3));
           gPos(i,2)=strtod(line(4));
           gPos(i,4)=strtod(line(10));
           if (line(5)=='+') then               
               gPos(i,3)=1;
           else
               gPos(i,3)=0;
           end
           if(int(strtod(line(10)))<690) then
               break
               end
           //break
       end
    end
    mclose(f);
end

fileName =proteinName+'_protein locations.txt';
f=mopen(fileName,'w');
for i=1:numOfgenes
    mputl((string(i)+ascii(9)+string(gPos(i,1))+ascii(9)+string(gPos(i,2))+ascii(9)+string(gPos(i,3))+ascii(9)+string(gPos(i,4))),f);
end
mclose(f);

fileName = proteinName+'_protein locations.txt';
gPos=zeros(numOfgenes,3);
f=mopen(fileName,'r');
for i=1:numOfgenes
    line=mgetl(f,1);
    indx=strindex(line,ascii(9));                            
    line=strsplit(line,indx);
    line=stripblanks(line,%t);
    gPos(i,1)=strtod(line(2));
    gPos(i,2)=strtod(line(3));
    gPos(i,3)=strtod(line(4));
end   
mclose(f)

FastaFolderPath = "C:\Users\Tharindu Rasanga\Documents\Maya\fasta"
f_files =  findfiles(FastaFolderPath,'*.fasta');
f_fileNames = gsort(f_files,'g','i');

geneSet=[];

genes=['NZ_CP014692.1','NC_013209.1/AP011121.1','NZ_CP023657.1/CP023657.1','NZ_CP022699.1/CP022699.1','NZ_CP014687.1/CP014687.1','NZ_LN606600.1/LN606600.1','NZ_CP011120.1/CP011120.1','NZ_CP015164.1','NZ_CP015168.1','NZ_CP021524.1','NZ_CP022374.1/CP022374.1','NZ_AP018515.1/AP018515.1','NZ_CP023189.1/CP023189.1','NC_017100.1/AP011128.1','NC_017121.1/AP011135.1','NC_017125.1/AP011142.1','NC_017146.1/AP011149.1','NZ_LN609302.1','NC_017111.1/AP011156.1','NC_017150.1/AP011163.1','NC_017108.1/AP011170.1','NZ_AP014881.1/AP014881.1']

for i=1:numOfgenes
    fileName=FastaFolderPath+"\"+f_fileNames(i);
    f=mopen(fileName,'r');
    disp(fileName);
    gene=get_fasta_at(fileName,gPos(i,1),gPos(i,2)+2,gPos(i,3));
    disp(ascii(gene));
    geneSet=[geneSet,ascii(gene)];
end

f=mopen('Set of extracted genes.txt','w');

for i=1:numOfgenes
    header='>'+genes(i);
    //header = ">" + ;
    mputl(header,f);
    mputl(geneSet(i),f);
end
mclose(f);

